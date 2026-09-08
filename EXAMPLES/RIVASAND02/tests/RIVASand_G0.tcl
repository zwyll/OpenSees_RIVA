# Both generations: required G0, units, element copies, elastic scale,
# nonlinear checkpoint continuation, and independence of material instances.
proc closeVector {a b label {rtol 1.0e-10}} {
    if {[llength $a] != [llength $b]} { error "$label: wrong length" }
    foreach x $a y $b {
        if {abs($x-$y)>$rtol*(1.0+abs($y))} { error "$label: $x != $y" }
    }
}
proc configure {} {
    wipeAnalysis
    constraints Transformation
    numberer RCM
    system BandGeneral
    test NormDispIncr 1.0e-12 30 0
    algorithm Newton
    integrator DisplacementControl 5 1 1.0e-7
    analysis Static
}
proc advance {first last} {
    set angle [expr {2.0*acos(-1.0)/32.0}]
    for {set i $first} {$i <= $last} {incr i} {
        set dg [expr {0.0002*(sin($i*$angle)-sin(($i-1)*$angle))}]
        integrator DisplacementControl 5 1 $dg
        if {[analyze 1] != 0} { error "G0 test failed at $i" }
    }
}
set g0ref 483.48301127222084
foreach name {RIVASand RIVASAND02} {
    foreach bad {0 -1 Inf NaN} {
        wipe
        model BasicBuilder -ndm 3 -ndf 3
        if {![catch {nDMaterial $name 1 .663 $bad 1.25 1.125 122.44207260468994 .945 .025 .78 .51 10 1.5 .65}]} {
            error "$name accepted invalid G0=$bad"
        }
    }
    if {![catch {nDMaterial $name 1 .663 1.25 1.125 122.44207260468994 .945 .025 .78 .51 10 1.5 .65}]} {
        error "$name accepted obsolete eleven-value input"
    }
    foreach g0 [list 350.0 $g0ref 700.0] {
        foreach scale {1.0 1000.0} {
            wipe
            model BasicBuilder -ndm 3 -ndf 3
            set initial [list [expr {-19.4*$scale}] [expr {-19.4*$scale}] [expr {-40.0*$scale}] 0 0 0]
            nDMaterial $name 1 .66296296296296309 $g0 1.25 1.125 122.44207260468994 .945 .025 .78 .51 10 1.5 .65 \
                -nSub 4 -stressScale $scale -stage 0 -initialStress {*}$initial
            # A second tag must not overwrite the first tag's stiffness.
            nDMaterial $name 2 .663 250.0 1.25 1.125 122.44207260468994 .945 .025 .78 .51 10 1.5 .65 -stressScale $scale
            foreach {tag x y z} {1 0 0 0 2 1 0 0 3 1 1 0 4 0 1 0 5 0 0 1 6 1 0 1 7 1 1 1 8 0 1 1} {
                node $tag $x $y $z
            }
            element bbarBrick 1 1 2 3 4 5 6 7 8 1
            foreach tag {1 2 3 4} { fix $tag 1 1 1 }
            foreach tag {5 6 7 8} { fix $tag 0 1 1 }
            foreach tag {6 7 8} { equalDOF 5 $tag 1 }
            pattern Plain 1 Linear { load 5 1.0 0.0 0.0 }
            closeVector [eleResponse 1 material 1 G0] $g0 "copied G0"
            set gref [expr {$g0*101.3*$scale}]
            closeVector [eleResponse 1 material 1 Gref] $gref "dimensional Gref"
            configure
            if {[analyze 1] != 0} { error "elastic test failed" }
            set stress [eleResponse 1 material 1 stress]
            closeVector [lindex $stress 5] [expr {$gref*1.e-7}] "elastic shear scale"
            updateMaterialStage -material 1 -stage 1
            advance 1 32
            set prefix [file join [pwd] G0_${name}_${g0}_${scale}]
            database File $prefix
            save 1
            advance 33 64
            set expected [eleResponse 1 material 1 stress]
            set expectedState [eleResponse 1 material 1 state]
            restore 1
            closeVector [eleResponse 1 material 1 G0] $g0 "restored G0"
            closeVector [eleResponse 1 material 1 Gref] $gref "restored Gref"
            configure
            advance 33 64
            closeVector [eleResponse 1 material 1 stress] $expected "restart stress" 1.e-12
            closeVector [eleResponse 1 material 1 state] $expectedState "restart state" 1.e-12
            if {$scale==1.0} { set kpaStress $expected } else {
                set normalized {}
                foreach x $expected { lappend normalized [expr {$x/$scale}] }
                closeVector $normalized $kpaStress "Pa/kPa equivalence" 1.e-8
            }
            wipe
            foreach path [glob -nocomplain ${prefix}*] { file delete -force $path }
        }
    }
}
foreach {emax emin} {.78 .51 .819 .51 .78 .49 .819 .49} {
    wipe
    model BasicBuilder -ndm 3 -ndf 3
    set drref [expr {($emax-.601)/($emax-$emin)}]
    nDMaterial RIVASAND02 1 $drref 700 1.25 1.125 122.44207260468994 .945 .025 $emax $emin 10 1.5 .65
    foreach {tag x y z} {1 0 0 0 2 1 0 0 3 1 1 0 4 0 1 0 5 0 0 1 6 1 0 1 7 1 1 1 8 0 1 1} {
        node $tag $x $y $z
        fix $tag 1 1 1
    }
    element bbarBrick 1 1 2 3 4 5 6 7 8 1
    closeVector [eleResponse 1 material 1 referenceRelativeDensity] $drref "custom reference density"
    closeVector [eleResponse 1 material 1 voidRatio] .601 "reference void-ratio anchor"
    set prefix [file join [pwd] G0_density_${emax}_${emin}]
    database File $prefix
    save 1
    restore 1
    closeVector [eleResponse 1 material 1 referenceRelativeDensity] $drref "restored reference density"
    closeVector [eleResponse 1 material 1 G0] 700 "restored custom density G0"
    wipe
    foreach path [glob -nocomplain ${prefix}*] { file delete -force $path }
}
puts "PASS: both RIVA-Sand generations, G0=350/reference/700, Pa/kPa, elastic scale, isolation, nonlinear restart, invalid and old inputs"
puts "PASS: material-specific reference density and G0 survive element copies/restart for four bound combinations"
