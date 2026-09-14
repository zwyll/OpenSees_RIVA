# Rotate a smooth combined strain path about the vertical axis at fixed rows.
load [lindex $argv 0]
set out [lindex $argv 1]
file mkdir $out
foreach {layer row} {
    bottom {0.90 900 1.25 1.125 122.44207260468994 0.945 0.025 0.78 0.51 4 4 0.65}
    loose {0.36 450 1.25 1.125 25 0.945 0.025 0.78 0.51 10 1.5 0.65}
} {
    foreach angle {0.0 0.61} {
        wipe
        model BasicBuilder -ndm 3 -ndf 4
        nDMaterial RIVASAND02BranchReversalResearch 1 {*}$row -stage 1 -initialStress -100 -100 -100 0 0 0 -nSub 4 -noBiasVolume
        foreach {n x y z} {1 0 0 0 2 1 0 0 3 1 1 0 4 0 1 0 5 0 0 1 6 1 0 1 7 1 1 1 8 0 1 1} {
            node $n $x $y $z
        }
        element SSPbrickUP 1 1 2 3 4 5 6 7 8 1 550000 1 1e-5 1e-5 1e-5 0.537 2.381e-6 0 0 0
        columnProbe update
        columnProbe commit
        set f [open "$out/${layer}_${angle}.csv" w]
        puts $f "time_s,sxx,syy,szz,sxy,syz,sxz,reversals"
        for {set i 1} {$i<=4096} {incr i} {
            set t [expr {$i/2048.0}]
            set gamma [expr {0.003*sin(2*acos(-1)*$t)}]
            foreach n {5 6 7 8} {
                setNodeDisp $n 1 [expr {$gamma*cos($angle)}] -commit
                setNodeDisp $n 2 [expr {$gamma*sin($angle)}] -commit
                setNodeDisp $n 3 [expr {0.0001*$t}] -commit
            }
            columnProbe update
            columnProbe commit
            puts $f "$t,[join [eleResponse 1 stress] ,],[eleResponse 1 reversals]"
        }
        close $f
    }
}
puts "ROTATION_COVARIANCE_COMPLETE"
wipe
