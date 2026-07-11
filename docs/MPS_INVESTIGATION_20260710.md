# Apple Metal GPU (MPS / MLX) usability from the sandboxed macOS compute environment

**Date:** 2026-07-10
**Investigator:** neoantigen-specialist subagent
**Question:** Why is the Apple Metal GPU not usable from this sandboxed macOS compute
environment, is it fixable from inside, and if not, what must change and who does it?

## Verdict (one line)

**Platform boundary.** The per-process Seatbelt sandbox denies the IOKit GPU user-client
open that Metal requires (`IOServiceOpen` → `kIOReturnNotPermitted`), so **zero** Metal
devices are exposed to any process in this environment. No env var, torch build, or MLX
version can fix it from inside; a platform/host change (launch the compute environment with
GPU/IOKit access, or run GPU work on Modal) is required.

## Environment (verified, not assumed)

| Fact | Value | How verified |
|---|---|---|
| Hardware | Apple M5 Max, 40 GPU cores, Metal: Supported | `system_profiler SPDisplaysDataType` |
| OS | macOS 26.5 (build 25F71), arm64 | `sw_vers`; `platform.mac_ver()` inside process → `('26.5', …)` |
| Kernel | Darwin 25.5.0 xnu-12377 T6050 | `uname -a` |
| User | `alex` uid=501, in `admin`/`_developer` groups | `id` |
| Session | **Aqua GUI login**, `hasGraphicAccess=True`, `isRemote=False` | `launchctl managername`; `Security.SessionGetInfo` attrs=0x6030 |
| Sandbox | **`SANDBOX_RUNTIME=1`** (Seatbelt active) | `env` |

The Aqua/`isRemote=False`/`hasGraphicAccess=True` result is important: this is **not** a
headless SSH session. The GPU denial therefore comes from the per-process sandbox, not from
the login-session layer.

## Reproduction (both failures, exact text)

**torch-MPS** (env `neoantigen-specialist`, torch 2.12.1):
```
torch 2.12.1
mps.is_built True
mps.is_available False
ALLOC ERR: RuntimeError The MPS backend is supported on macOS 14.0+.
           Current OS version can be queried using `sw_vers`
```
The "macOS 14.0+" string is a **misleading fallback message**, not a real version gate: the
process itself reads the OS as 26.5 (`platform.mac_ver()`). torch prints this whenever it
cannot obtain a Metal device, regardless of OS version.

**MLX** (`pip install mlx`, current):
```
RuntimeError: [metal::load_device] No Metal device available. This typically occurs in
headless, sandboxed, or virtualized macOS sessions where the GPU is not accessible.
```
MLX's own diagnostic names the cause: a sandboxed session where the GPU is not accessible.

## Root cause (evidence chain)

The failure is **below** the torch/MLX layer. Four independent probes, none using
torch or MLX, all show zero Metal devices:

1. **Python ctypes → Metal framework directly**
   - `MTLCreateSystemDefaultDevice()` → **NULL**
   - `MTLCopyAllDevices()` → array returned but **count = 0**
   - The Metal framework loads and the calls succeed; the process is simply exposed no GPUs.

2. **Native compiled C** (`clang -framework Metal`)
   - `objc default=nil copyall_count=0` — confirms it is not a Python/ctypes artifact.

3. **Distinction from the ordinary headless quirk.** In a plain headless-SSH macOS session,
   `MTLCreateSystemDefaultDevice()` returns nil but `MTLCopyAllDevices()` **still enumerates**
   the GPU. Here **both** report zero. This is a stronger isolation than headless: the device
   is not merely "not default", it is not enumerable at all.

4. **IOKit user-client probe (the mechanism).** The GPU driver *is* visible in the IO
   registry — `IOServiceGetMatchingServices` matches `IOGPU`, `AGXAccelerator`, and
   `IOAccelerator`, each returning a found service. But opening the user-client, which is what
   Metal must do internally to create a device, fails:
   ```
   IOGPU          match_kr=0 first_service=found  IOServiceOpen_kr=-536870206
   AGXAccelerator match_kr=0 first_service=found  IOServiceOpen_kr=-536870206
   IOAccelerator  match_kr=0 first_service=found  IOServiceOpen_kr=-536870206
   ```
   `-536870206 = 0xE00002C2 = kIOReturnNotPermitted`.

**Conclusion:** hardware and driver are present and visible to the process; the Seatbelt
sandbox **denies the IOKit user-client open** (`kIOReturnNotPermitted`) that Metal requires
to instantiate a device. Corroborating sandbox evidence: `swiftc` could not write its module
cache ("Operation not permitted"), and there are no GPU `/dev` nodes. This is category **(b):
a platform/sandbox boundary where the Metal device is not exposed to this process**, not
category (a) a fixable software/config problem.

## Candidate fixes tested

| # | Candidate fix | Result | Interpretation |
|---|---|---|---|
| 1 | `PYTORCH_ENABLE_MPS_FALLBACK=1` | `default=NULL copyall_count=0` | No effect (only affects op fallback *after* a device exists) |
| 2 | `MTL_DEVICE_WRAPPER_TYPE=1` | `default=NULL copyall_count=0` | No effect |
| 3 | `METAL_DEVICE_WRAPPER_TYPE=1` | `default=NULL copyall_count=0` | No effect |
| 4 | `MTL_HUD_ENABLED=1` | `default=NULL copyall_count=0` | No effect |
| 5 | `MTL_DEBUG_LAYER=0 MTL_SHADER_VALIDATION=0` | `default=NULL copyall_count=0` | No effect |
| 6 | `CA_ALLOW_HEADLESS=1` | `default=NULL copyall_count=0` | No effect |
| 7 | `__GL_SYNC_TO_VBLANK=0` | `default=NULL copyall_count=0` | No effect |
| 8 | Baseline (no vars) | `default=NULL copyall_count=0` | Reference |
| 9 | Alternate torch build / nightly | **Not applicable** | Native compiled C (no torch) also gets nil; the boundary is below the torch layer, so no torch build can create a device the OS refuses to hand this process |
| 10 | MLX version bump | **Not applicable** | Same reason — MLX calls the same `MTLCreateSystemDefaultDevice`; current MLX already reproduces the failure |
| 11 | Newer OS (torch "14.0+" message) | **Not applicable** | Process already sees macOS 26.5; the message is a fallback string, not a real gate |

No in-sandbox configuration changed the outcome. Every path bottoms out at the same
`kIOReturnNotPermitted` IOKit denial.

## What would NOT be a legitimate fix (explicitly not attempted)

The sandbox GPU isolation is a security boundary. The following were **not** attempted and
should not be: disabling Seatbelt, granting the process a GPU IOKit entitlement by
re-signing, sandbox escape, running the interpreter outside the sandbox profile, or any
privilege escalation. Any of these would be circumvention, not a fix.

## Remediation path

This is **not fixable from inside the sandbox**. To get Apple-GPU compute, one of:

1. **Platform/host action (owner of this compute environment):** launch the macOS compute
   sandbox with a profile that permits the GPU IOKit user-client — i.e. allow `IOServiceOpen`
   on `IOGPU`/`AGXAccelerator` so `MTLCreateSystemDefaultDevice()` succeeds. This is a
   host-launch / sandbox-profile change made by the platform, **not** a code change the agent
   or user can make from within the session. Until that is done, MPS and MLX will remain
   unavailable here regardless of code.

2. **Use the configured remote GPU fallback (recommended, available now):** run Evo 2 /
   HyenaDNA and any MPS/MLX/CUDA work on **Modal (`byoc:modal`)**, which provides real GPU
   passthrough. No platform change to this Mac sandbox is required.

**Recommendation for the project:** treat local Apple-GPU compute as unavailable in this
environment and route all genomic-foundation-model scoring to Modal. Local CPU remains fine
for lightweight work (a `pytorch cpuonly` env, `seqenc`, already exists).

## Working recipe (in-sandbox)

None. There is no in-sandbox configuration that yields a usable Metal device.
