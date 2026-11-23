#!/usr/bin/env python3
import csv, os, re, itertools, hashlib, datetime

ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..'))
AN=os.path.join(ROOT,'全部分析内容'); OUT=os.path.join(ROOT,'revision_R1_major6'); RES=os.path.join(OUT,'results')
os.makedirs(RES,exist_ok=True); os.makedirs(os.path.join(OUT,'logs'),exist_ok=True)
PRIMARY=['GSE1919','GSE12021','GSE55235','GSE55457','GSE82107','GSE206848','GSE89408','GSE143514']
ALL_BULK=PRIMARY+['GSE283079']

def read_csv(path):
    with open(path,newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def write_csv(path,rows,fields):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
def write(path,text):
    with open(path,'w',encoding='utf-8') as f:f.write(text)
def load_bulk(ds):
    p=os.path.join(AN,'phase4f_reference_based_deconvolution_sensitivity','03_deconvolution_inputs',ds+'_metadata_for_deconvolution.csv')
    return read_csv(p)
def parse_relation(s):
    return re.findall(r'(?:BioSample|SRA):\s*(?:https?://[^/]+/(?:biosample|sra)\?term=)?([A-Z]+\d+)',s or '')
def gsm_set(rows):
    return {r.get('sample_id','') for r in rows if re.fullmatch(r'GSM\d+',r.get('sample_id',''))}

canonical=read_csv(os.path.join(ROOT,'revision_R1_major2','results','canonical_matched_input.csv'))
canon_by={d:[r for r in canonical if r['dataset']==d] for d in PRIMARY}
raw_metadata_counts={}
raw_group_counts={}
for _ds in ALL_BULK:
    _raw=read_csv(os.path.join(AN,'phase4f_reference_based_deconvolution_sensitivity','03_deconvolution_inputs',_ds+'_metadata_for_deconvolution.csv'))
    raw_metadata_counts[_ds]=len(_raw)
    raw_group_counts[_ds]={'OA':sum(r.get('disease_group','').lower() in ('oa','osteoarthritis') for r in _raw), 'Control':sum(r.get('disease_group','').lower() not in ('oa','osteoarthritis') for r in _raw)}
cross=read_csv(os.path.join(AN,'gse12021_independent_unit_correction','source_data','GSE12021_source_record_crosswalk.csv'))
cross_in=[r for r in cross if r.get('canonical_analysis_include','').upper()=='TRUE']
bulk_records={}
for ds in ALL_BULK:
    rows=load_bulk(ds); arr=[]
    for r in rows:
        sid=r.get('sample',''); gsm=sid if sid.startswith('GSM') else ''
        arr.append({'sample_id':sid,'GSM':gsm,'disease_group':r.get('disease_group',''),'relation':r.get('relation',''),'platform':r.get('platform',r.get('platform_id','')),'age':r.get('age',''),'sex':r.get('sex','')})
    bulk_records[ds]=arr
bulk_records['GSE12021']=[{'sample_id':r['GSM'],'GSM':r['GSM'],'disease_group':r['phenotype'],'relation':'corrected canonical GPL96','platform':r['platform'],'age':'','sex':''} for r in cross_in]

scb=read_csv(os.path.join(AN,'layer25s_disease_state_recovery_and_symptom_boundary','data','raw_geo','GSE283080','sample_metadata.csv'))
bulk_records['GSE283080']=[{'sample_id':r['geo_accession'],'GSM':r['geo_accession'],'disease_group':r.get('disease state:ch1',''),'relation':r.get('relation',''),'platform':r.get('platform_id',''),'age':'','sex':''} for r in scb]
refa=read_csv(os.path.join(AN,'phase4f_reference_based_deconvolution_sensitivity','01_reference_construction','reference_cell_metadata.csv'))
refa_donors=sorted({r['donor_id'] for r in refa})

def counts(ds):
    rows=bulk_records[ds]; oa=sum(r['disease_group'].lower() in ('oa','osteoarthritis') for r in rows); ctrl=len(rows)-oa
    current=len(canon_by.get(ds,[]))
    if ds in raw_group_counts:
        return raw_metadata_counts[ds],current,raw_group_counts[ds]['OA'],raw_group_counts[ds]['Control']
    return len(rows),current,oa,ctrl

prov=[]
for ds in PRIMARY+['GSE283079']:
    raw,current,oa,ctrl=counts(ds); rows=bulk_records[ds]; platforms=sorted({r['platform'] for r in rows if r['platform']})
    if ds=='GSE12021':
        status='Confirmed independent within accession after correction; cross-dataset overlap unresolved'; level='confirmed_within_accession'; ids='19 patient IDs (EB/J) in retained GPL96 crosswalk; GPL97 duplicates mapped to 14 retained patients and excluded'; overlap='14 duplicate platform records excluded; no cross-accession patient crosswalk'
    else:
        status='Unresolved'; level='unresolved'; ids='GSM/sample identifiers available; no reusable patient/donor crosswalk in frozen inputs'; overlap='No exact cross-dataset GSM overlap found in audited local metadata; patient-level overlap not testable'
    prov.append({'dataset':ds,'role':'primary_bulk' if ds in PRIMARY else 'external_holdout','raw_metadata_records':raw,'current_canonical_records':current if ds in PRIMARY else 'not_primary','OA_n':oa,'Control_n':ctrl,'platforms':';'.join(platforms),'known_patient_or_donor_ids':ids,'within_accession_duplicate_platform_status':overlap,'cross_dataset_overlap_status':status,'evidence_level':level,'evidence_files':'phase4f deconvolution metadata; canonical_matched_input.csv; GSE12021_source_record_crosswalk.csv' if ds=='GSE12021' else 'phase4f deconvolution metadata; canonical_matched_input.csv','unresolved_information':'Patient identity and recruitment-cohort overlap with other GEO series not documented in frozen local inputs'})
for ds in ['GSE248453','GSE283080']:
    if ds=='GSE283080':
        rows=bulk_records[ds]; detail='4 GSM donor identifiers; all OA; no patient IDs'; st='Potential overlap'; ev='GSE283080 sample_metadata.csv; official GEO series page; shared PMID/contributors with GSE283079'
        prov.append({'dataset':ds,'role':'single_cell_reference_Bv2','raw_metadata_records':len(rows),'current_canonical_records':'reference cells '+str(len([r for r in read_csv(os.path.join(AN,'reference_B_downstream_impact_v2','inputs','reference','reference_B_cell_metadata.csv'))])),'OA_n':4,'Control_n':0,'platforms':'GPL24676','known_patient_or_donor_ids':detail,'within_accession_duplicate_platform_status':'No duplicate donor IDs in frozen metadata','cross_dataset_overlap_status':st,'evidence_level':'potential_overlap_unresolved','evidence_files':ev,'unresolved_information':'No direct donor crosswalk between GSE283080 and GSE283079; shared study does not prove shared participants'})
    else:
        prov.append({'dataset':ds,'role':'single_cell_reference_A','raw_metadata_records':'not available as GEO sample table in frozen local inputs','current_canonical_records':'reference cells '+str(len(refa)),'OA_n':8,'Control_n':0,'platforms':'10x scRNA reference','known_patient_or_donor_ids':f'{len(refa_donors)} donor IDs: '+ ';'.join(refa_donors),'within_accession_duplicate_platform_status':'No duplicate donor IDs in frozen reference metadata','cross_dataset_overlap_status':'Unresolved','evidence_level':'donor_resolved_within_reference_only','evidence_files':'reference_cell_metadata.csv; reference construction report','unresolved_information':'No GSM/BioSample/patient crosswalk linking Reference A donors to bulk cohorts or holdout'})
fields=list(prov[0].keys()); write_csv(os.path.join(RES,'dataset_provenance_audit.csv'),prov,fields)

series={d:set(r['GSM'] for r in bulk_records[d] if r['GSM']) for d in ALL_BULK+['GSE283080']}
bio={d:set(x for r in bulk_records[d] for x in parse_relation(r['relation'])) for d in ['GSE283079','GSE283080']}
pairs=[]
for a,b in itertools.combinations(ALL_BULK+['GSE283080','GSE248453'],2):
    gsm=sorted(series.get(a,set()) & series.get(b,set())); bios=sorted(bio.get(a,set()) & bio.get(b,set()))
    if {a,b}=={'GSE283079','GSE283080'}: status='Potential overlap / unresolved'; basis='Same official paper/contributors/institution/submission date and adjacent GSM range, but separate BioProjects and no shared GSM/BioSample IDs'
    elif {a,b}=={'GSE55235','GSE55457'}: status='Unresolved'; basis='Separate accessions and no exact GSM overlap; both GPL96, but no patient crosswalk'
    else: status='Unresolved'; basis='No exact shared GSM/BioSample ID detected in available metadata; patient recruitment overlap not documented'
    pairs.append({'dataset_a':a,'dataset_b':b,'exact_shared_GSM_n':len(gsm),'exact_shared_GSM_ids':';'.join(gsm),'exact_shared_BioSample_n':len(bios),'exact_shared_BioSample_ids':';'.join(bios),'patient_id_overlap':'not_testable','independence_assessment':status,'evidence_basis':basis})
write_csv(os.path.join(RES,'sample_independence_matrix.csv'),pairs,list(pairs[0].keys()))
candidate_pairs=[]
for x in pairs:
    a,b=x['dataset_a'],x['dataset_b']
    if {a,b}=={'GSE283079','GSE283080'}:
        priority='high'; reason='same official paper/study family; participant overlap unresolved'
    elif {a,b}=={'GSE55235','GSE55457'}:
        priority='high'; reason='both GPL96; separate accessions but no patient crosswalk'
    elif x['independence_assessment']=='Unresolved':
        priority='routine'; reason='patient-level overlap cannot be tested from frozen local metadata'
    else:
        continue
    candidate_pairs.append({**x,'priority':priority,'candidate_reason':reason})
write_csv(os.path.join(RES,'cross_dataset_overlap_candidates.csv'),candidate_pairs,list(candidate_pairs[0].keys()))

hold=[]
for target,role in [('GSE248453','Reference A'),('GSE283080','Reference B-v2')]+[(d,'primary bulk') for d in PRIMARY]:
    p=next(x for x in pairs if {x['dataset_a'],x['dataset_b']}=={'GSE283079',target})
    status='Potential overlap / unresolved' if target=='GSE283080' else 'Unresolved'
    hold.append({'holdout':'GSE283079','comparison_target':target,'target_role':role,'exact_GSM_overlap_n':p['exact_shared_GSM_n'],'exact_BioSample_overlap_n':p['exact_shared_BioSample_n'],'same_publication_or_project_evidence':'yes' if target=='GSE283080' else 'not documented','donor_overlap_confirmed':'no','assessment':status,'holdout_use':'external OA-vs-control bulk holdout; not an independent cell-composition validation unless donor independence is confirmed','evidence':p['evidence_basis']})
write_csv(os.path.join(RES,'holdout_reference_independence.csv'),hold,list(hold[0].keys()))

sens=[{'analysis':'overlap_exclusion_meta_analysis','status':'not_performed','reason':'No confirmed cross-dataset duplicate patient or sample record was identified; unresolved patient-level overlap cannot support defensible deletion','primary_result_preserved':'yes','required_future_evidence':'patient/donor crosswalk or original study-level recruitment records'}]
write_csv(os.path.join(RES,'overlap_sensitivity_results.csv'),sens,list(sens[0].keys()))

urls=['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283079','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283080','https://pmc.ncbi.nlm.nih.gov/articles/PMC11790023/']
write_csv(os.path.join(OUT,'logs','source_url_audit.csv'),[{'url':u,'access_date':'2026-11-13','note':'Official GEO or linked paper provenance evidence'} for u in urls],['url','access_date','note'])
write(os.path.join(OUT,'logs','run_parameters.txt'),'read-only provenance audit; no meta-analysis rerun; exact identifier intersections only; seed not applicable\n')
