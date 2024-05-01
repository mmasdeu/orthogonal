$(document).ready(function () {
	//Only needed for the filename of export files.
	//Normally set in the title tag of your page.
	document.title = "DGL Table";
	// Create search inputs in footer
	$("#dgltable tfoot th").each(function () {
		var title = $(this).text();
		$(this).html('<input type="text" placeholder="Search ' + title + '" />');
	});
	// DataTable initialisation
	var table = $("#dgltable").DataTable({
		dom: '<"dt-buttons"Bf><"clear">lirtp',
		paging: false,
		autoWidth: true,
		buttons: [
			"colvis",
			"copyHtml5",
			"csvHtml5",
		],
		initComplete: function (settings, json) {
			var footer = $("#dgltable tfoot tr");
			$("#dgltable thead").append(footer);
		}
	});

	// Apply the search
	$("#dgltable thead").on("keyup", "input", function () {
		table.column($(this).parent().index())
		.search(this.value,true,false)
		.draw();
	});
});
