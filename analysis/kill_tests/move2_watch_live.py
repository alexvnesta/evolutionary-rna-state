"""Live in-session watch: poll the store for a changed quant_gene_tpm.parquet and
auto-fire move2_autorun. Bounded so it self-terminates. Prints one line per poll
so liveness is visible; writes results as artifacts are handled by the caller."""
import importlib.util, json, time, sys

def watch(host, poll_s=120, max_hours=6):
    spec = importlib.util.spec_from_file_location("m2", "move2_autorun.py")
    M2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(M2)
    deadline = time.time() + max_hours * 3600
    n_fired = 0
    last = None
    while time.time() < deadline:
        try:
            r = M2.run_if_changed(host)
        except Exception as e:
            print(f"[watch] poll error: {type(e).__name__}: {e}", flush=True)
            time.sleep(poll_s); continue
        st = r.get("status")
        fp = (r.get("version_id"), r.get("n_samples"))
        if st == "ran":
            n_fired += 1
            print(f"[watch] FIRED on {r.get('n_samples')} samples "
                  f"{r.get('cohorts')} -> {r.get('out_prefix')} "
                  f"(T2 p={r.get('T2_p')}, activity AUROC={r.get('T3_activity_auroc')})",
                  flush=True)
        elif fp != last:
            print(f"[watch] {st}: {r.get('msg','')} ({r.get('n_samples')} samples)", flush=True)
        last = fp
        time.sleep(poll_s)
    return {"fired": n_fired, "stopped": "deadline"}

if __name__ == "__main__":
    try:
        host  # noqa: F821
    except NameError:
        print("run inside kernel", file=sys.stderr); sys.exit(2)
    print(json.dumps(watch(host)))
