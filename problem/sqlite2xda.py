import sqlite3
from typing import IO, Optional

NOSCALE = ("year", "stand", )

def outcda(db: sqlite3.Connection, out: IO, sample: Optional[int] = None):
	standquery = f"""
		SELECT COUNT(*), stand
		FROM tail
		WHERE stand%? = 0
		GROUP BY stand
		ORDER BY stand
	"""
	for (ns, id) in db.execute(standquery, (sample or 1, )):
		out.write(f"{ns} {int(id)}\n")

def outxda(db: sqlite3.Connection, xda: IO, cols: list[str], sample: Optional[int] = None):
	tailname = [v for v in cols if "#" not in v]
	tailidx = [i for i,v in enumerate(cols) if "#" not in v]
	nodename = list({v.split("#")[0] for v in cols if "#" in v})
	periods = max(int(v.split("#")[1]) for v in cols if "#" in v)
	nodeidx = [[
		next((i for i,v in enumerate(cols) if v==f"{n}#{p+1}"), None)
		for n in nodename
	] for p in range(periods)]
	tailquery = f"""
		SELECT id {''.join(f',{v}' for v in tailname)}
		FROM tail
		WHERE stand%? = 0
		ORDER BY stand
	"""
	nodequery = f"""
		WITH branch AS (
			SELECT parent {''.join(f',{v}' for v in nodename)}
				FROM nodes
				WHERE id=?
			UNION ALL
			SELECT n.parent {''.join(f',n.{v}' for v in nodename)}
				FROM branch b, nodes n
				WHERE n.id=b.parent
		)
		SELECT {','.join(nodename)}
			FROM branch
			ORDER BY parent
	"""
	scale = [c not in NOSCALE and sample or 1 for c in cols]
	for (id, *tail) in db.execute(tailquery, (sample or 1, )):
		xdarow = [None] * len(cols)
		for i,v in zip(tailidx, tail):
			xdarow[i] = v
		prev = [0] * len(nodename)
		for (p, node) in enumerate(db.execute(nodequery, (id, ))):
			for i,v,p in zip(nodeidx[p], node, prev):
				if i is not None:
					xdarow[i] = v-p
			prev = node
		xda.write(" ".join(f"{x*c:.10g}" for x,c in zip(xdarow, scale))) # type: ignore
		xda.write("\n")

if __name__ == "__main__":
	import os
	import sys
	try:
		in_ = sys.argv[1]
		out = sys.argv[2]
		if not (out.endswith(".cda") or out.endswith(".xda")):
			raise Exception
	except:
		print(f"usage:\n  {sys.argv[0]} data.sqlite3 out.cda\n  {sys.argv[0]} data.sqlite3 out.xda var1 ... varN")
		sys.exit(1)	
	scale = int(os.environ["XDA_SCALE"]) if "XDA_SCALE" in os.environ else None
	db = sqlite3.connect(in_)
	with open(out, "w") as fp:
		if out.endswith(".cda"):
			outcda(db, fp, scale)
		else:
			vars = [v.strip() for v in open("../variables.txt")]
			outxda(db, fp, vars, scale)
