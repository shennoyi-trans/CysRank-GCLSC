param([switch]$CheckOnly, [switch]$NoBrowser, [int]$Port = 8765)
$ErrorActionPreference = 'Stop'
[Console]::OutputEncoding = New-Object System.Text.UTF8Encoding($false)
$OutputEncoding = [Console]::OutputEncoding
$root = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $root
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
# Ignore inherited Python overrides; always use a verified executable explicitly.
Remove-Item Env:PYTHONHOME -ErrorAction SilentlyContinue
Remove-Item Env:PYTHONPATH -ErrorAction SilentlyContinue
$env:PYTHONNOUSERSITE = '1'
$logDir = Join-Path $root 'logs\launcher'
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$transcript = Join-Path $logDir (Get-Date -Format 'yyyyMMdd-HHmmss-fff').ToString()
Start-Transcript -Path "$transcript.log" | Out-Null
$startupLock = $null
try {
    # Hold an exclusive project lock through setup and serving: a second double-click
    # must not run concurrent pip installs against the same environment.
    try {
        $startupLock = [IO.File]::Open((Join-Path $logDir 'startup.lock'), 'OpenOrCreate', 'ReadWrite', 'None')
    } catch { throw 'Another launcher is running for this project. Use its browser window, or close it before restarting.' }
    if (-not [Environment]::Is64BitOperatingSystem -or $env:PROCESSOR_ARCHITECTURE -eq 'ARM64' -or $env:PROCESSOR_ARCHITEW6432 -eq 'ARM64') {
        throw 'This launcher requires Windows x64 (Intel/AMD).'
    }
    Write-Host 'Cyclic Intelligence: checking Python and installed dependencies...'
    $candidates = @(
        (Join-Path $root '.venv\Scripts\python.exe'),
        (Join-Path $root 'tmp\archive\.venv-gpu\Scripts\python.exe'),
        (Join-Path $root '.runtime\python313\python.exe')
    )
    $candidates += @(Get-ChildItem -Path (Join-Path $root '.runtime\venv-*\Scripts\python.exe') -ErrorAction SilentlyContinue | ForEach-Object { $_.FullName })
    $candidates += @(Get-Command python.exe -All -ErrorAction SilentlyContinue | ForEach-Object { $_.Source })
    if (Get-Command py.exe -ErrorAction SilentlyContinue) {
        $candidates += @(& py.exe -0p 2>$null | ForEach-Object {
            if ($_ -match '([A-Za-z]:\\.*python(?:3\.13)?\.exe)\s*$') { $Matches[1] }
        })
    }
    foreach ($key in @('HKCU:\Software\Python\PythonCore\3.13\InstallPath', 'HKLM:\Software\Python\PythonCore\3.13\InstallPath')) {
        if (Test-Path $key) { $candidates += (Join-Path (Get-Item $key).GetValue('') 'python.exe') }
    }
    $condaList = Join-Path $env:USERPROFILE '.conda\environments.txt'
    if (Test-Path -LiteralPath $condaList) {
        $candidates += @(Get-Content -LiteralPath $condaList | ForEach-Object { Join-Path $_ 'python.exe' })
    }
    $best = $null
    $bestScore = -1
    foreach ($candidate in ($candidates | Select-Object -Unique)) {
        if (-not (Test-Path -LiteralPath $candidate) -or $candidate -like '*\Microsoft\WindowsApps\*') { continue }
        try {
            $probeText = & $candidate (Join-Path $PSScriptRoot 'launch.py') --probe 2>$null
            if ($LASTEXITCODE -ne 0) { continue }
            $probe = $probeText | ConvertFrom-Json
            if ($probe.compatible -and $probe.score -gt $bestScore) {
                $best = $candidate
                $bestScore = $probe.score
            }
        } catch { continue }
    }
    if (-not $best) {
        if ($CheckOnly) { throw 'No usable CPython 3.13 x64 found. Run start.cmd without -CheckOnly to install it.' }
        $runtime = Join-Path $root '.runtime'
        New-Item -ItemType Directory -Force -Path $runtime | Out-Null
        $installer = Join-Path $runtime 'python-3.13.5-amd64.exe'
        Write-Host 'Downloading official Python 3.13.5 x64 installer...'
        [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
        if (-not (Test-Path -LiteralPath $installer)) {
            Invoke-WebRequest -UseBasicParsing -Uri 'https://www.python.org/ftp/python/3.13.5/python-3.13.5-amd64.exe' -OutFile "$installer.part"
            Move-Item -LiteralPath "$installer.part" -Destination $installer
        }
        $signature = Get-AuthenticodeSignature -LiteralPath $installer
        if ($signature.Status -ne 'Valid' -or $signature.SignerCertificate.Subject -notmatch 'Python Software Foundation') {
            throw "Python installer signature verification failed: $installer"
        }
        $target = Join-Path $runtime 'python313'
        Write-Host 'Installing Python for the current user (no PATH changes)...'
        $installArgs = "/quiet InstallAllUsers=0 TargetDir=`"$target`" Include_pip=1 Include_launcher=0 Include_test=0 Include_doc=0 Include_tcltk=0 Shortcuts=0 PrependPath=0 AssociateFiles=0"
        $process = Start-Process -FilePath $installer -ArgumentList $installArgs -WindowStyle Hidden -Wait -PassThru
        if ($process.ExitCode -notin @(0, 3010)) { throw "Python installer failed with exit code $($process.ExitCode)." }
        $best = Join-Path $target 'python.exe'
        if (-not (Test-Path -LiteralPath $best)) { throw 'Python installation did not produce the expected interpreter.' }
    }
    Write-Host "Selected Python: $best"
    $launchArgs = @('-u', (Join-Path $PSScriptRoot 'launch.py'), '--port', "$Port")
    if ($CheckOnly) { $launchArgs += '--check-only' }
    if ($NoBrowser) { $launchArgs += '--no-browser' }
    # PowerShell 5.1 transcripts omit native output unless it is sent through the host.
    # Native stderr is logged, but only the process exit code determines failure.
    $savedErrorAction = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try {
        & $best @launchArgs 2>&1 | ForEach-Object { Write-Host $_ }
        $launchExitCode = $LASTEXITCODE
    } finally { $ErrorActionPreference = $savedErrorAction }
    if ($launchExitCode -ne 0) { throw "Launcher exited with code $launchExitCode." }
} catch {
    Write-Host "`nERROR: $($_.Exception.Message)" -ForegroundColor Red
    Write-Host "Log: $transcript.log"
    exit 1
} finally {
    if ($startupLock) { $startupLock.Dispose() }
    Stop-Transcript | Out-Null
}
