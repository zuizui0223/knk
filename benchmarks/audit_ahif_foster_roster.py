"""Outcome-blind source taxon roster bridge for published 41-species sowing list."""
import argparse,csv,json
from pathlib import Path
from benchmarks.audit_ahif_seed_synthesis import acquire
from benchmarks.score_ahif_seed_synthesis_opposition import read
from benchmarks.score_ahif_kansas_forecast import EXPERIMENT,grid

def norm(s):
    return " ".join(s.split()).casefold()

def audit_names(names,rows):
    p={norm(x["paper_species"]):x for x in rows}
    q={norm(n):n for n in names}
    if len(p)!=41 or len(q)!=100 or len(rows)!=41 or len(names)!=100:raise ValueError("name roster collision or count")
    exact=[{"source":q[n],"paper":p[n]["paper_species"],"guild":p[n]["functional_guild"],
            "reported_present_before_sowing":p[n]["reported_present_before_sowing"]} for n in sorted(p.keys()&q.keys())]
    return {"exact_matches":exact,"unmatched_source_names":[q[n] for n in sorted(q.keys()-p.keys())],
            "unmatched_paper_names":[p[n]["paper_species"] for n in sorted(p.keys()-q.keys())],
            "n_exact":len(exact),"n_paper":len(p),"n_source":len(q)}

def main():
    a=argparse.ArgumentParser();a.add_argument("--root",type=Path,required=True);a.add_argument("--out",type=Path,required=True);a.add_argument("--roster",type=Path,required=True);args=a.parse_args()
    args.root.mkdir(parents=True,exist_ok=True)
    m=acquire(args.root);s,roster,dup,drift=read(args.root,m)
    if dup[EXPERIMENT] or drift[EXPERIMENT]:raise ValueError("physical key ambiguity")
    x,blocks=grid(s[EXPERIMENT])
    with args.roster.open(newline="",encoding="utf-8") as f:rows=list(csv.DictReader(f))
    if sum(x["reported_present_before_sowing"]=="yes" for x in rows)!=3:raise ValueError("preexisting roster")
    z={"status":"PUBLISHED_ROSTER_BRIDGE_NAMES_ONLY_NOT_ORIGINAL_PLOT_AUDIT",
       "paper_doi":"10.1890/08-0849.1","source_doi":"10.6084/m9.figshare.12319682.v1",
       "source_original_cited_doi":"10.6073/pasta/5703e8819c81da0beed4660f467c288e",
       "n_states":120,"n_blocks":len(blocks),"n_ages":15,**audit_names(roster[EXPERIMENT],rows),
       "limits":["Probable shared source experiment, original EDI subset not independently verified.",
                 "No synonym mapping or response calculation yet. Membership is assigned seed list, not individual origin."]}
    args.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print("FOSTER_ROSTER_AUDIT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
