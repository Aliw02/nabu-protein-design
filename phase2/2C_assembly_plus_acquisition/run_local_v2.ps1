$ErrorActionPreference = "Stop"

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\\..")
Set-Location $repoRoot

$dataDir = Join-Path $repoRoot ".phase2c_data"
New-Item -ItemType Directory -Force -Path $dataDir | Out-Null

$source = Join-Path $dataDir "gb1_wu2016_landscape.csv"
$url = "https://raw.githubusercontent.com/rcronk/mutate/0861a3f1f6c58a9597fde689eb7f9a07b38e2638/ridge/data/gb1_wu2016_landscape.csv"
$expectedSha = "02ef9bc1784a67ee14cf4024596a8fcfc7129779911c6589ed896373552125b1"

Write-Host "[2C] Downloading pinned GB1 source..."
Invoke-WebRequest -Uri $url -OutFile $source

$actualSha = (Get-FileHash -Path $source -Algorithm SHA256).Hash.ToLowerInvariant()
if ($actualSha -ne $expectedSha) {
    throw "GB1 SHA256 mismatch. Expected $expectedSha, got $actualSha"
}

Write-Host "[2C] Source hash verified."
Write-Host "[2C] Running preregistered V2 full-beam maturity experiment..."

$outDir = "phase2/2C_assembly_plus_acquisition/results/gb1_gated_v2_local"
python phase2/2C_assembly_plus_acquisition/run_gb1_2c.py $source --probe-count 256 --out $outDir

if ($LASTEXITCODE -ne 0) {
    throw "Phase 2C V2 runner failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "[2C] Comparison:"
Get-Content (Join-Path $repoRoot "$outDir/PHASE2C_COMPARISON.csv")

Write-Host ""
Write-Host "[2C] Manifest:"
Write-Host (Join-Path $repoRoot "$outDir/PHASE2C_MANIFEST.json")
