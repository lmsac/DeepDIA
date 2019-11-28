Param($fastaPath, $accessionPath)

function Filter-FastaByAccession() {

    param(
        [string]$fastaPath,
        [string]$accessionPath
    )

    [string]$ext = [IO.Path]::GetExtension($fastaPath)  
    [string]$outPath = [IO.Path]::ChangeExtension($fastaPath, "filtered$ext")

    $accessionSet = New-Object System.Collections.Generic.HashSet[string]

    Get-Content $accessionPath | ForEach-Object {
        $accession = $_.Trim()
        if (!$accessionSet.Contains($accession)) {
            $accessionSet.Add($accession)
        }
    } | Out-Null

    [IO.StreamReader]$reader = [IO.File]::OpenText($fastaPath)
    [System.IO.StreamWriter]$writer = [IO.File]::CreateText($outPath)
     
    try {
                
        [Text.StringBuilder]$title = New-Object Text.StringBuilder   
        [Text.StringBuilder]$sequence = New-Object Text.StringBuilder
        while($true) {   
            if (!$reader.EndOfStream) {         
                $line = $reader.ReadLine()
            }
            if ($reader.EndOfStream -or $line.StartsWith('>')) {                
                if ($title.Length -gt 0) {  
                    $accession = $title.ToString().Substring(1)
                    if (!$accessionSet.Contains($accession)) {
                        $accession = $accession.Split(' ')[0].Trim()
                    }
                    if (!$accessionSet.Contains($accession)) {
                        $accession = $accession.Split('|')[1]
                    }
                    if ($accessionSet.Contains($accession)) {
			$writer.Write($title)
                        $writer.Write($sequence)
			# Write-Host "$title"
                    }
                }

                $title.Clear().AppendLine($line) | Out-Null
                $sequence.Clear() | Out-Null
            }
            else {
                $sequence.AppendLine($line) | Out-Null
            }  
            if ($reader.EndOfStream) {
                break
            }          
        }
    }
    finally {
        $reader.Close()
        $writer.Close()
    }
}

Filter-FastaByAccession $fastaPath $accessionPath