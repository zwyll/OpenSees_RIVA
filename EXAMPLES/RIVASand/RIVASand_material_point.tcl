# RIVA-Sand zero-bias reference, exercised through a homogeneous 3D brick.
# Stress unit: kPa. The imposed top displacement equals engineering gamma_xz
# because the specimen height is one.

wipe
model BasicBuilder -ndm 3 -ndf 3

set matTag 8001
set Dr [expr {(0.78-0.601)/(0.78-0.51)}]
nDMaterial RIVASand $matTag \
    $Dr 1.25 1.125 122.44207260468994 0.945 0.025 \
    0.78 0.51 10.0 1.5 0.65 \
    -nSub 1 -stressScale 1.0 \
    -initialStress -19.4 -19.4 -40.0 0.0 0.0 0.0

# Unit cube, with node ordering used by the OpenSees brick elements.
node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0
element bbarBrick 1 1 2 3 4 5 6 7 8 $matTag

# Enforce homogeneous simple shear: the bottom is fixed and all top nodes
# share one horizontal displacement.  Vertical and transverse movements are
# held to zero so the prescribed strain path matches the material-point oracle.
fix 1 1 1 1
fix 2 1 1 1
fix 3 1 1 1
fix 4 1 1 1
fix 5 0 1 1
fix 6 0 1 1
fix 7 0 1 1
fix 8 0 1 1
equalDOF 5 6 1
equalDOF 5 7 1
equalDOF 5 8 1

pattern Plain 1 Linear {
    load 5 1.0 0.0 0.0
}

constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-12 12 0
algorithm Newton
integrator DisplacementControl 5 1 1.0e-12
analysis Static

set output [open "RIVASand_material_point.csv" "w"]
puts $output "step,cycle,gamma_xz,tau_xz,sigma_xx,sigma_yy,sigma_zz,skeleton_pressure_loss_ratio,reversals,compatibility_residual"

set pi [expr {acos(-1.0)}]
set pointsPerCycle 32
set totalSteps [expr {2*$pointsPerCycle}]
set amplitude 0.005
set previousGamma 0.0

for {set step 0} {$step <= $totalSteps} {incr step} {
    set cycle [expr {double($step)/$pointsPerCycle}]
    set gamma [expr {$amplitude*sin(2.0*$pi*$cycle)}]
    if {$step > 0} {
        set increment [expr {$gamma-$previousGamma}]
        integrator DisplacementControl 5 1 $increment
        if {[analyze 1] != 0} {
            error "RIVASand brick test failed at step $step"
        }
    }
    set stress [eleResponse 1 material 1 stress]
    set pressureLoss [lindex [eleResponse 1 material 1 effectivePressureRatio] 0]
    set reversals [lindex [eleResponse 1 material 1 reversals] 0]
    set residual [lindex [eleResponse 1 material 1 compatibilityResidual] 0]
    puts $output [join [list $step $cycle $gamma \
        [lindex $stress 5] [lindex $stress 0] [lindex $stress 1] \
        [lindex $stress 2] $pressureLoss $reversals $residual] ","]
    set previousGamma $gamma
}

close $output
puts "Wrote RIVASand_material_point.csv"
