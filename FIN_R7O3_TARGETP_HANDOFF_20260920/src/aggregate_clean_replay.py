from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'checkpoints/clean_math_shards';TOTAL=12425
rows=[];cover=[];failures=[]
for p in sorted(D.glob('shard_*.json')):
 d=json.load(open(p));rows.append({'file':p.name,**{k:d[k] for k in ['range','count','pass','fail','seconds']}});cover.extend(range(d['range'][0],d['range'][1]));failures+=d.get('failures',[])
missing=sorted(set(range(TOTAL))-set(cover));dups=len(cover)-len(set(cover));extra=sorted(set(cover)-set(range(TOTAL)))
out={'task':'R7O3-clean-directory-full-math-replay','shard_count':len(rows),'covered_count':len(set(cover)),'pass_count':sum(r['pass'] for r in rows),'fail_count':sum(r['fail'] for r in rows),'duplicate_index_count':dups,'missing':missing,'extra':extra,'failure_examples':failures[:50],'global_pass':len(set(cover))==TOTAL and not missing and not extra and dups==0 and sum(r['pass'] for r in rows)==TOTAL and sum(r['fail'] for r in rows)==0,'shards':rows}
(ROOT/'results').mkdir(exist_ok=True);(ROOT/'results/R7O3-037_clean_directory_math_replay.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['shards','failure_examples']},indent=2));sys.exit(0 if out['global_pass'] else 1)
