template_text = """$$comment
$comment$$end

$$molecule
$charge $multiplicity
$geometry$$end

$$rem
jobtype $job_type
exchange $exchange
basis $basis
molden_format $molden_format
$other_rems
$$end

$other_things"""

# other things is for stuff like $basis or $ecp etc..
