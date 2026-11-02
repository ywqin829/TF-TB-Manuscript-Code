$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$OutputPath = Join-Path $ProjectRoot "docs\DELIVERY_MANIFEST.tsv"
$ExcludedPrefix = Join-Path $ProjectRoot "data\GSE192483_3-Celltype"

$Rows = Get-ChildItem -LiteralPath $ProjectRoot -Recurse -File |
    Where-Object {
        $_.FullName -ne $OutputPath -and
        -not $_.FullName.StartsWith(
            $ExcludedPrefix,
            [System.StringComparison]::OrdinalIgnoreCase
        )
    } |
    ForEach-Object {
        $rootWithSeparator = $ProjectRoot.TrimEnd("\") + "\"
        $relative = $_.FullName.Substring($rootWithSeparator.Length)
        $top = ($relative -split '[\\/]')[0]
        $role = switch ($top) {
            "submission_package" { "Final upload-ready submission package" }
            "scripts" { "Reproducible analysis or document-build code" }
            "results" { "Analysis result, cache, or intermediate table" }
            "figures" { "Main, supplementary, or source-context figure" }
            "docs" { "Manuscript, cover letter, review, guide, or QC report" }
            "supplement" { "Submission supplementary tables or index" }
            "logs" { "Audit log, session record, or internal traceability" }
            default { "Project documentation" }
        }
        [pscustomobject]@{
            RelativePath = $relative
            AbsolutePath = $_.FullName
            Role = $role
            SizeBytes = $_.Length
            LastWriteTime = $_.LastWriteTime.ToString("s")
        }
    } |
    Sort-Object RelativePath

$Rows | Export-Csv -LiteralPath $OutputPath -Delimiter "`t" -NoTypeInformation -Encoding UTF8
