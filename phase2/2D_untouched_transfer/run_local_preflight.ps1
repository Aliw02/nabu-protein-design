$ErrorActionPreference = "Stop"

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\\..")
Set-Location $repoRoot

$dataDir = Join-Path $repoRoot ".phase2d_blind_source"
$outDir = Join-Path $repoRoot ".phase2d_local"
New-Item -ItemType Directory -Force -Path $dataDir | Out-Null
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

$aavZip = Join-Path $dataDir "aav_splits.zip"
$aavReference = Join-Path $dataDir "P03135.fasta"
$ired = Join-Path $dataDir "ired_two_to_many.csv.gz"

$aavUrl = "https://raw.githubusercontent.com/microsoft/protein-uq/5e7b2b9cd219805eaabe73b21d4d9955cf882448/data/aav/splits.zip"
$aavReferenceUrl = "https://raw.githubusercontent.com/microsoft/protein-uq/5e7b2b9cd219805eaabe73b21d4d9955cf882448/data/aav/P03135.fasta"
$iredUrl = "https://flip.protein.properties/assets/splits/ired/two_to_many.csv.gz"

$expectedAav = "ad91ba8d5b390d793fc9393f8003ff7b0290fbe2e13db58cb6ad72bc981a99fd"
$expectedAavReference = "9b5e572c2a18d27482b629efeb48f573e866fe61bb6aea58bb8c087e4177fbcb"
$expectedIred = "aa45a2f85fb1af87b6b0e86397b3f8f292a061dc574be2536673b5bd490b0e74"

Write-Host "[2D] Running contract tests..."
python -m pytest -q phase2/2D_untouched_transfer/tests
if ($LASTEXITCODE -ne 0) { throw "Phase 2D contract tests failed." }

Write-Host "[2D] Downloading pinned sources..."
Invoke-WebRequest -Uri $aavUrl -OutFile $aavZip
Invoke-WebRequest -Uri $aavReferenceUrl -OutFile $aavReference
Invoke-WebRequest -Uri $iredUrl -OutFile $ired

function Assert-Sha256([string]$Path, [string]$Expected) {
    $actual = (Get-FileHash -Path $Path -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($actual -ne $Expected) { throw "SHA256 mismatch for $Path. Expected $Expected, got $actual" }
    Write-Host "[2D] SHA256 OK: $Path"
}

Assert-Sha256 $aavZip $expectedAav
Assert-Sha256 $aavReference $expectedAavReference
Assert-Sha256 $ired $expectedIred

Write-Host "[2D] Running identity-only preflight..."
python phase2/2D_untouched_transfer/preflight_identity.py --aav-zip $aavZip --aav-reference $aavReference --ired-gz $ired --out (Join-Path $outDir "IDENTITY_PREFLIGHT.json")
if ($LASTEXITCODE -ne 0) { throw "Phase 2D identity preflight failed." }

Write-Host "[2D] Running IRED scoreability preflight..."
python phase2/2D_untouched_transfer/ired_scoreability_preflight.py --ired-gz $ired --out (Join-Path $outDir "IRED_SCOREABILITY_PREFLIGHT.json")
if ($LASTEXITCODE -ne 0) { throw "IRED scoreability preflight failed." }

Write-Host "[2D] LABEL-BLIND PREFLIGHT COMPLETE"
Write-Host "[2D] Identity report: $outDir\\IDENTITY_PREFLIGHT.json"
Write-Host "[2D] IRED scoreability report: $outDir\\IRED_SCOREABILITY_PREFLIGHT.json"