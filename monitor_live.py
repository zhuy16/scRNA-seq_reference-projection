#!/usr/bin/env python3
"""
Live resource monitor — run this in a separate terminal while training.

Usage:  python3 monitor_live.py
        python3 monitor_live.py --pid <PID>   # attach to specific process
        python3 monitor_live.py --interval 2  # sample every 2 s (default: 3)

Press Ctrl-C to stop.
"""
import argparse, os, sys, time, psutil

try:
    import torch
    MPS = torch.backends.mps.is_available()
except ImportError:
    torch = None
    MPS = False

def mps_mem():
    if MPS:
        try:
            return torch.mps.current_allocated_memory() / 1073741824
        except Exception:
            pass
    return 0.0

def find_training_pids():
    """Find Python processes that might be doing SCVI training."""
    pids = []
    for p in psutil.process_iter(['pid', 'name', 'cmdline']):
        try:
            if 'python' in p.info['name'].lower():
                cmd = ' '.join(p.info['cmdline'] or [])
                if 'ipykernel' in cmd or 'scvi' in cmd.lower() or 'papermill' in cmd:
                    pids.append(p.info['pid'])
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return pids

def fmt_bar(pct, width=20):
    filled = int(pct / 100 * width)
    return '[' + '#' * filled + '.' * (width - filled) + ']'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pid', type=int, default=None)
    ap.add_argument('--interval', type=float, default=3.0)
    args = ap.parse_args()

    print(f"\n{'='*62}")
    print(f"  Live Resource Monitor  (Ctrl-C to stop, interval={args.interval}s)")
    print(f"  MPS (Apple GPU): {'available' if MPS else 'not available'}")
    print(f"{'='*62}")

    proc = None
    if args.pid:
        proc = psutil.Process(args.pid)
        print(f"  Monitoring PID {args.pid}: {' '.join(proc.cmdline()[:3])}\n")

    header = (f"  {'Time':>8}  {'Sys CPU%':>9}  {'SysFree GB':>11}"
              f"  {'ProcRAM GB':>11}  {'MPS GB':>7}")
    sep    = "  " + "-"*8 + "  " + "-"*9 + "  " + "-"*11 + "  " + "-"*11 + "  " + "-"*7

    print(header)
    print(sep)

    t0 = time.time()
    peak_cpu = 0.0
    peak_ram = 0.0
    min_free = float('inf')
    peak_mps = 0.0

    try:
        while True:
            vm      = psutil.virtual_memory()
            cpu     = psutil.cpu_percent(interval=None)
            free_gb = vm.available / 1073741824
            mps_gb  = mps_mem()
            elapsed = time.time() - t0

            proc_ram = 0.0
            if proc:
                try:
                    proc_ram = proc.memory_info().rss / 1073741824
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            else:
                # sum all python/ipykernel processes
                for p in find_training_pids():
                    try:
                        proc_ram += psutil.Process(p).memory_info().rss / 1073741824
                    except Exception:
                        pass

            peak_cpu = max(peak_cpu, cpu)
            peak_ram = max(peak_ram, proc_ram)
            min_free = min(min_free, free_gb)
            peak_mps = max(peak_mps, mps_gb)

            bar_cpu  = fmt_bar(cpu)
            bar_free = fmt_bar(100 - vm.percent, 10)

            ts = f"{int(elapsed//60):02d}:{int(elapsed%60):02d}"
            print(f"  {ts:>8}  {cpu:8.1f}%  {bar_free} {free_gb:4.2f}GB"
                  f"  {proc_ram:10.2f}GB  {mps_gb:6.3f}GB")
            sys.stdout.flush()
            time.sleep(args.interval)

    except KeyboardInterrupt:
        print(f"\n{sep}")
        print(f"  {'PEAK':>8}  {peak_cpu:8.1f}%  {'':11}  {peak_ram:10.2f}GB  {peak_mps:6.3f}GB")
        print(f"  Min system free RAM: {min_free:.2f} GB")
        print(f"  Total monitored time: {int((time.time()-t0)//60)}m {int((time.time()-t0)%60)}s")
        print(f"{'='*62}\n")

if __name__ == '__main__':
    main()
