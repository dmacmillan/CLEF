$(document).ready(function() {

// Loading samples 
// **************************************************************
$('#load_sample').click(function() {
    if ($('#runbutton').length) {
	    patients = '001	A02:01:01G	A03:01:01G	B0801G	B51:01:01G	C04:01:01G	C14:02:01G	ATGTGCTGCGGGATCCGAGAC\n' + 
        '002	A02:01:01G	A02:01:01G	B35:02:01G	B51:01:01G	C04:01:01G	C16:02:01G	ATGCGCTGCGGCTTCCGAGAC';
    } else {
	    patients = '001	A02:01:01G	A03:01:01G	B35:01:01G	B51:01:01G	C04:01:01G	C14:02:01G\n' + 
        '002	A02:01:01G	A02:01:01G	B35:02:01G	B51:01:01G	C04:01:01G	C16:02:01G';
        sequences = 'ATGTGCTGCGGGATCCGAGAC\n' + 'ATGCGCTGCGGCTTCCGAGAC'
	    $('#sequences').val(sequences);
	    $('#sequences').text(sequences);
    }
   patients = '001	A02:01:01G	A03:01:01G	B08:01:01G	B51:01:01G	C04:01:01G	C14:02:01G	ATGTGCTGCGGGATCCGAGAC\n' + 
        '002	A02:01:01G	A02:01:01G	B01:01:01G	B51:01:01G	C04:01:01G	C16:02:01G	ATGCGCTGCGGCTTCCGAGAC';
	hlas = 'B*08:01	2R	nonadapted\nB*08:01	2C	adapted'
	$('#hlas').text(hlas);
	//$('#patients').text(patients);
	$('#hlas').val(hlas);
   $('#patients').text(patients)
   $('#patients').val(patients)
	//$('#patients').val(patients);
});
// **************************************************************

// Clearing
// **************************************************************
$('#clear').click(function() {
	$('#patients').text('');
	$('#patients').val('');
	$('#patients').focus();
});
// **************************************************************

$('.clearbutton').click(function() {
	$(this).siblings('textarea').text('');
	$(this).siblings('textarea').val('');
});

});
