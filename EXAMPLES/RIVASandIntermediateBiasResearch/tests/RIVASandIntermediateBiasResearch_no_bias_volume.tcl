# Verify the opt-in no-bias-volume adapter flag and restart continuity.

proc makeBrick {matTag} {
    foreach {tag x y z} {
        1 0 0 0  2 1 0 0  3 1 1 0  4 0 1 0
        5 0 0 1  6 1 0 1  7 1 1 1  8 0 1 1
    } {
        node $tag $x $y $z
        fix $tag 1 1 1
    }
    element bbarBrick 1 1 2 3 4 5 6 7 8 $matTag
}

set Dr [expr {(0.78-0.601)/(0.78-0.51)}]
set databasePrefix [file join [pwd] RIVASandIntermediateBiasResearch_no_bias_db]
foreach path [glob -nocomplain ${databasePrefix}*] { file delete -force $path }

# The accepted research response remains the default.
wipe
model BasicBuilder -ndm 3 -ndf 3
nDMaterial RIVASandIntermediateBiasResearch 8020 \
    $Dr 1.25 1.125 122.44207260468994 0.945 0.025 \
    0.78 0.51 10.0 1.5 0.65 -stressScale 1.0
makeBrick 8020
set defaultFlag [lindex [eleResponse 1 material 1 noBiasVolume] 0]
if {$defaultFlag != 0.0} {
    error "RIVASandIntermediateBiasResearch unexpectedly disables bias volume by default"
}

set conflictRejected [catch {
    nDMaterial RIVASandIntermediateBiasResearch 8022 \
        $Dr 1.25 1.125 122.44207260468994 0.945 0.025 \
        0.78 0.51 10.0 1.5 0.65 -stressScale 1.0 \
        -fieldBiasVolume -noBiasVolume
}]
if {!$conflictRejected} {
    error "-fieldBiasVolume and -noBiasVolume were incorrectly accepted together"
}

# The option must survive the element's material copy and database restart.
wipe
model BasicBuilder -ndm 3 -ndf 3
nDMaterial RIVASandIntermediateBiasResearch 8021 \
    $Dr 1.25 1.125 122.44207260468994 0.945 0.025 \
    0.78 0.51 10.0 1.5 0.65 -stressScale 1.0 -noBiasVolume
makeBrick 8021
set copiedFlag [lindex [eleResponse 1 material 1 noBiasVolume] 0]
if {$copiedFlag != 1.0} {
    error "element material copy lost -noBiasVolume"
}

database File $databasePrefix
save 201
restore 201
set restoredFlag [lindex [eleResponse 1 material 1 noBiasVolume] 0]
if {$restoredFlag != 1.0} {
    error "sendSelf/recvSelf lost -noBiasVolume"
}

wipe
foreach path [glob -nocomplain ${databasePrefix}*] { file delete -force $path }
puts "PASS: -noBiasVolume is opt-in and survives material copy and save/restore"
