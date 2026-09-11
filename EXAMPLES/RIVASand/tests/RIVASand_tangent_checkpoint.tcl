# Run with the pre-change executable in "make" mode, then the new executable
# in "check" mode, in the SAME scratch directory. Tests revision-1 -> 2 reads.
set ::rivaTangentLibraryOnly 1
source [file join [file dirname [info script]] RIVASand_continuum.tcl]
if {[llength $argv]!=1 || [lindex $argv 0] ni {make check}} {error "expected make or check"}
foreach name {RIVASand RIVASAND02} {
    set prefix [file join [pwd] old_tangent_${name}]
    set expectedFile ${prefix}.expected
    if {[lindex $argv 0] eq "make"} {
        setup $name -1 1 .663 .25 4 0
        advance 20 0
        database File $prefix
        save 1
        set expected [lindex [advance 44 0] 0]
        set fp [open $expectedFile w]; puts $fp $expected; close $fp
    } else {
        wipe; model BasicBuilder -ndm 3 -ndf 3
        database File $prefix
        restore 1
        closeVector [eleResponse 1 material 1 TanType] 0 "old checkpoint defaults to elastic" 0
        configure
        set actual [lindex [advance 44] 0]
        set fp [open $expectedFile r]; set expected [read $fp]; close $fp
        closeVector [concat {*}$actual] [concat {*}$expected] "old checkpoint continuation" 0
    }
}
puts "PASS: checkpoint [lindex $argv 0]"
