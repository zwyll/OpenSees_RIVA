# Prescribed smooth shear reversals, with and without simultaneous vertical strain.
load [lindex $argv 0]
set out [lindex $argv 1]
file mkdir $out
set summary [open "$out/summary.csv" w]
puts $summary "material,vertical_rate,steps_per_cycle,reversals,final_tau_kPa"
foreach material {RIVASAND02 RIVASAND02BranchReversalResearch} {
    foreach vertical {0.0 0.0001} {
        foreach steps {128 512 2048 8192} {
            wipe
            model BasicBuilder -ndm 3 -ndf 4
            nDMaterial $material 1 0.90 900 1.25 1.125 122.44207260468994 0.945 0.025 0.78 0.51 4 4 0.65 -stage 1 -initialStress -100 -100 -100 0 0 0 -nSub 4 -noBiasVolume
            foreach {n x y z} {1 0 0 0 2 1 0 0 3 1 1 0 4 0 1 0 5 0 0 1 6 1 0 1 7 1 1 1 8 0 1 1} {
                node $n $x $y $z
            }
            element SSPbrickUP 1 1 2 3 4 5 6 7 8 1 550000 1 1e-5 1e-5 1e-5 0.537 2.381e-6 0 0 0
            columnProbe update
            columnProbe commit
            set f [open "$out/${material}_v${vertical}_n${steps}.csv" w]
            puts $f "time_s,gamma_xz,epsilon_zz,tau_xz_kPa,pressure_kPa,reversals"
            for {set i 1} {$i<=[expr {2*$steps}]} {incr i} {
                set t [expr {double($i)/$steps}]
                set gamma [expr {0.003*sin(2*acos(-1)*$t)}]
                foreach n {5 6 7 8} {
                    # Commit node coordinates only so the second assignment
                    # retains the first component; material commit follows update.
                    setNodeDisp $n 1 $gamma -commit
                    setNodeDisp $n 3 [expr {$vertical*$t}] -commit
                }
                columnProbe update
                columnProbe commit
                set stress [eleResponse 1 stress]
                set pressure [expr {-([lindex $stress 0]+[lindex $stress 1]+[lindex $stress 2])/3}]
                puts $f "$t,$gamma,[expr {$vertical*$t}],[lindex $stress 5],$pressure,[eleResponse 1 reversals]"
            }
            close $f
            puts $summary "$material,$vertical,$steps,[eleResponse 1 reversals],[lindex [eleResponse 1 stress] 5]"
            flush $summary
        }
    }
}
close $summary
puts "ROTATING_PATH_CHECKS_COMPLETE"
wipe
