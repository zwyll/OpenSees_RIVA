# Small boundary-value convergence probe, NOT a LEAP/full-mesh qualification.
# Usage: OpenSees .../RIVASand_tangent_benchmark.tcl
# Same stress update, fixed dt, nSub, loads, tolerance and solver in each pair.
foreach name {RIVASand RIVASAND02} {
    foreach kind {staticRamp staticLow staticHigh dynamicUP biasedUP} {
        foreach type {0 1} {
            wipe
            set up [expr {$kind in {dynamicUP biasedUP}}]
            set bias [expr {$kind eq "biasedUP"?10.0:0.0}]
            model BasicBuilder -ndm 3 -ndf [expr {$up?4:3}]
            set opts {}
            if {$name eq "RIVASAND02"} { lappend opts -reversalLatch }
            nDMaterial $name 1 .663 483.48301127222084 1.25 1.125 122.44207260468994 .945 .025 .78 .51 10 1.5 .65 \
                -TanType $type -nSub 4 -rho 1.8 -stage 0 -initialStress -19.4 -19.4 -40 0 0 $bias {*}$opts
            foreach {tag x y z} {1 0 0 0 2 1 0 0 3 1 1 0 4 0 1 0 5 0 0 1 6 1 0 1 7 1 1 1 8 0 1 1} {
                node $tag $x $y $z
                if {$z==0} {
                    if {$up} {fix $tag 1 1 1 0} else {fix $tag 1 1 1}
                } else {
                    if {$up} {fix $tag 0 1 1 0} else {fix $tag 0 1 1}
                }
            }
            foreach tag {6 7 8} {equalDOF 5 $tag 1}
            if {$up} {
                # Undrained, finite fluid compressibility; no fluid stiffness
                # belongs in the NDMaterial tangent. No Rayleigh damping.
                element brickUP 1 1 2 3 4 5 6 7 8 1 2.2e6 1.0 1.e-8 1.e-8 1.e-8 0 0 0
            } else {element bbarBrick 1 1 2 3 4 5 6 7 8 1}
            updateMaterialStage -material 1 -stage 1
            constraints Transformation
            numberer RCM
            system BandGeneral
            test NormDispIncr 1.e-10 50 0
            algorithm Newton
            set values {0}
            set dt [expr {$up?.01:1.0}]
            set steps [expr {$up?300:128}]
            set amplitude [expr {$kind eq "staticLow"?1.0:4.0}]
            for {set i 1} {$i<=$steps} {incr i} {
                set cycles [expr {$up?$i*.01:$i/32.}]
                lappend values [expr {$kind eq "staticRamp"?$amplitude*$i/$steps:$amplitude*sin(2*acos(-1)*$cycles)}]
            }
            timeSeries Path 1 -dt $dt -values {*}$values -useLast
            if {$up} {pattern Plain 1 1 {load 5 1 0 0 0}} else {pattern Plain 1 1 {load 5 1 0 0}}
            if {$bias!=0} {pattern Plain 2 Constant {load 5 $bias 0 0 0}}
            if {$up} {integrator Newmark .5 .25; analysis Transient} else {integrator LoadControl $dt; analysis Static}
            set count 0; set iterations 0; set active 0; set code 0
            set started [clock microseconds]
            for {set i 0} {$i<$steps} {incr i} {
                if {$up} {set code [analyze 1 $dt]} else {set code [analyze 1]}
                incr iterations [testIter]
                if {$code!=0} {break}
                incr count
                if {[lindex [eleResponse 1 material 1 tangentStatus] 0]==1} {incr active}
            }
            set seconds [expr {([clock microseconds]-$started)/1.e6}]
            puts "BENCH $name $kind TanType=$type completed=$count/$steps iterations=$iterations continuumSteps=$active seconds=$seconds code=$code"
        }
    }
}
