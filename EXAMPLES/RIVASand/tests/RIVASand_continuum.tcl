# Both adapters: prescribed affine strains, tangent queries, units, copying,
# zero-increment trials, all BiasVolume modes, and database continuation.
# Run in a scratch directory; an optional output filename records the mode-0
# fingerprint so the pre-change executable can be compared independently.
proc closeVector {a b label {tol 1.e-11}} {
    if {[llength $a] != [llength $b]} { error "$label: length" }
    foreach x $a y $b {
        if {abs($x-$y)>$tol*(1+abs($y))} { error "$label: $x != $y" }
    }
}
proc setup {name type scale dr bias nsub mode {floor .5065}} {
    wipe
    model BasicBuilder -ndm 3 -ndf 3
    set options {}
    if {$type>=0} { lappend options -TanType $type }
    if {$name eq "RIVASAND02"} { lappend options -BiasVolume $mode -reversalLatch }
    nDMaterial $name 1 $dr 483.48301127222084 1.25 1.125 122.44207260468994 .945 .025 .78 .51 10 1.5 .65 \
        -nSub $nsub -stressScale $scale -stage 1 -tangentPMin [expr {$floor*$scale}] \
        -initialStress [expr {-19.4*$scale}] [expr {-19.4*$scale}] [expr {-40*$scale}] 0 0 [expr {$bias*40*$scale}] {*}$options
    foreach {tag x y z} {1 0 0 0 2 1 0 0 3 1 1 0 4 0 1 0 5 0 0 1 6 1 0 1 7 1 1 1 8 0 1 1} {
        node $tag $x $y $z
        if {$z==0} { fix $tag 1 1 1 } else { fix $tag 0 1 1 }
    }
    element bbarBrick 1 1 2 3 4 5 6 7 8 1
    # A separate elastic DOF avoids a zero-equation system; the brick strain
    # is prescribed, not determined by the selected material tangent.
    node 101 2 0 0; node 102 3 0 0
    fix 101 1 1 1; fix 102 0 1 1
    uniaxialMaterial Elastic 2 100
    element truss 2 101 102 1 2
    set values {0.0}
    for {set i 1} {$i<=128} {incr i} {
        lappend values [expr {.003*sin(2*acos(-1)*$i/32.)}]
    }
    timeSeries Path 1 -dt 1 -values {*}$values
    pattern Plain 1 1 { foreach node {5 6 7 8} { sp $node 1 1 } }
    pattern Plain 2 Linear { load 102 1 0 0 }
    configure
}
proc configure {} {
    wipeAnalysis
    constraints Transformation
    numberer RCM
    system BandGeneral
    test NormDispIncr 1.e-12 20
    algorithm Newton
    integrator LoadControl 1
    analysis Static
}
proc advance {count {queries 1}} {
    set history {}; set active 0; set mappingActive 0
    for {set i 0} {$i<$count} {incr i} {
        if {[analyze 1]!=0} { error "prescribed-strain test failed" }
        set state [eleResponse 1 material 1 state]
        set stress [eleResponse 1 material 1 stress]
        if {$queries} {
            set tangent [eleResponse 1 material 1 tangent]
            if {[llength $tangent]!=36} { error "tangent matrix missing" }
            foreach v $tangent { if {abs($v)>1.e100} { error "invalid tangent" } }
            set status [lindex [eleResponse 1 material 1 tangentStatus] 0]
            if {$status==1} {
                incr active
                if {[llength $state]==139 && [lindex $state 125]>0} {incr mappingActive}
            }
            # Querying a tangent must never advance the constitutive state.
            closeVector [eleResponse 1 material 1 state] $state "query state" 0
        }
        lappend history [concat $stress $state]
    }
    return [list $history $active $mappingActive]
}
if {[info exists ::rivaTangentLibraryOnly]} {return}
if {[llength $argv]} {
    # Supported by the old executable: no new options or tangent responses.
    set fp [open [lindex $argv 0] w]
    foreach name {RIVASand RIVASAND02} {
        setup $name -1 1 .663 .25 4 0
        puts $fp [lindex [advance 128 0] 0]
    }
    close $fp
    puts "PASS: default fingerprint written"
    exit
}
foreach name {RIVASand RIVASAND02} {
    foreach tail {{-TanType -1} {-TanType 2} {-TanType .5} {-TanType NaN} {-TanType Inf} {-TanType} {-TanType 0 -TanType 1}} {
        wipe; model BasicBuilder -ndm 3 -ndf 3
        if {![catch {nDMaterial $name 1 .663 483.48301127222084 1.25 1.125 122.44207260468994 .945 .025 .78 .51 10 1.5 .65 {*}$tail}]} {
            error "$name accepted invalid option $tail"
        }
    }
    set modes {0}
    if {$name eq "RIVASAND02"} { set modes {0 1 2} }
    foreach mode $modes {
        foreach {dr bias} {.663 0.0 .663 .25 .5 .25} {
            foreach nsub {1 4 16} {
                foreach type {-1 0 1} {
                    setup $name $type 1 $dr $bias $nsub $mode
                    closeVector [eleResponse 1 material 1 TanType] [expr {max($type,0)}] "copied selector" 0
                    lassign [advance 64] history active mappingActive
                    if {$type==-1} { set baseline $history } else {
                        closeVector [concat {*}$history] [concat {*}$baseline] "stress/state unchanged" 0
                    }
                    if {$type==1 && $active<1} { error "$name $dr $bias nsub=$nsub continuum never active" }
                    if {$name eq "RIVASAND02" && $type==1 && $dr==.5 && $bias==.25 && $mappingActive<1} {
                        error "mapping-backstress continuum branch not exercised"
                    }
                }
            }
        }
    }
    foreach mode $modes { foreach scale {1 1000} {
        setup $name 1 $scale .663 .25 4 $mode
        advance 20
        set matrix [eleResponse 1 material 1 tangent]
        set status [eleResponse 1 material 1 tangentStatus]
        set prefix [file join [pwd] tangent_${name}_${mode}_${scale}]
        database File $prefix
        save 1
        set expected [lindex [advance 44] 0]
        restore 1
        closeVector [eleResponse 1 material 1 TanType] 1 "restored selector" 0
        if {$name eq "RIVASAND02"} {closeVector [eleResponse 1 material 1 BiasVolume] $mode "restored volume mode" 0}
        closeVector [eleResponse 1 material 1 tangentStatus] $status "restored status" 0
        closeVector [eleResponse 1 material 1 tangent] $matrix "restored tangent" 0
        configure
        closeVector [concat {*}[lindex [advance 44] 0]] [concat {*}$expected] "restart continuation" 0
        if {$scale==1} { set kpa $matrix } else {
            set normalized {}; foreach v $matrix { lappend normalized [expr {$v/$scale}] }
            closeVector $normalized $kpa "Pa/kPa continuum" 1.e-8
        }
        wipe
        foreach path [glob -nocomplain ${prefix}*] { file delete $path }
    } }
    setup $name 1 1 .663 0 4 0 100
    lassign [advance 32] history active
    if {$active!=0} { error "tangent-pressure floor safeguard missing" }
}
puts "PASS: continuum options, defaults, prescribed-strain parity, copies, queries, restart, units, mapping and floors"
