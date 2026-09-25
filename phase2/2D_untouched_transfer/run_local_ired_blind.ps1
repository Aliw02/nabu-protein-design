$ErrorActionPreference = "Stop"

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\\..")
Set-Location $repoRoot

$dataDir = Join-Path $repoRoot ".phase2d_blind_source"
$outDir = Join-Path $repoRoot ".phase2d_local\\ired_blind_v2"
$ired = Join-Path $dataDir "ired_two_to_many.csv.gz"
$expectedIred = "aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74"

if (-not (Test-Path $ired)) { throw "Pinned IRED source not found. Run run_local_preflight.ps1 first." }
$actual = (Get-FileHash -Path $ired -Algorithm SHA256).Hash.ToLowerInvariant()
if ($actual -ne $expectedIred) { throw "IRED SHA256 mismatch. Expected $expectedIred, got $actual" }

Write-Host "[2D] Re-running contract tests before target reveal..."
python -m pytest -q phase2/2D_untouched_transfer/tests
if ($LASTEXITCODE -ne 0) { throw "Phase 2D contract tests failed. Blind target reveal aborted." }

Write-Host "[2D] STARTING IRED BLIND TARGET REVEAL"
Write-Host "[2D] Frozen policy: IRED_ABSTENTION_V2"
python phase2/2D_untouched_transfer/run_ired_blind_v2.py $ired --out $outDir
if ($LASTEXITCODE -ne 0) { throw "IRED blind run failed. Preserve this failed run before any fix." }

Write-Host "[2D] IRED BLIND RUN COMPLETE"
Write-Host "[2D] Manifest: $outDir\\IRED_BLIND_MANIFEST.json"
Write-Host "[2D] Predictions: $outDir\\IRED_TEST_PREDICTIONS.csv"