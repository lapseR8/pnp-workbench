import subprocess, sys, os
shard, n = sys.argv[1], sys.argv[2]
subprocess.run("curl -s https://data.pequalnp.com/data/desk/gaps.json -o gaps.json", shell=True, check=True)
os.environ["SHARD"], os.environ["NSHARDS"] = shard, n
subprocess.run([sys.executable, "upbatch.py", "gaps.json", "4"])
