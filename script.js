
$(document).ready(function () {
	//Only needed for the filename of export files.
	//Normally set in the title tag of your page.
    document.title = "DGL Table";
    // DataTable initialisation
    $('#dgltable').on( 'dblclick', 'tbody td', function (e) {
	if (e.ctrlKey) {
    let data = table.row(e.target.closest('tr')).data();
	let sage_code = "";
	let prec = Math.ceil((Math.log(data[12]) / Math.log(data[0])));
	if (prec == 0) {
		prec = 100;
	}
	if (data[5] == 'triv') {
		sage_code = "prime_list = [" + data[9].replace(/<[^>]*>?/gm, '').split("[")[1].split("]")[0] + "]; p = " + data[0] + "; prec = " + prec + "; label = " + data[1] + "; D = " + data[2] + "; n = " + data[3] + "; type = '" + data[4] + "'; char = '" + data[5] + "'; J0 = " + data[12] + "; J = Qq(p**2, prec, names='t')(2*J0);";
	}
	else {
		sage_code = "prime_list = [" + data[9].replace(/<[^>]*>?/gm, '').split("[")[1].split("]")[0] + "]; p = " + data[0] + "; prec = " + prec + "; label = " + data[1] + "; D = " + data[2] + "; n = " + data[3] + "; type = '" + data[4] + "'; char = '" + data[5] + "'; J0 = " + data[12] + "; J = Qq(p**2, prec, names='t')(J0) + ((1-Qq(p**2, prec, names='t')(J0)**2)/(t - t.trace()/2).norm()).sqrt() * (t - t.trace()/2);";
	}
	
	navigator.clipboard.writeText(sage_code);

	// Alert the copied text
	alert("Copied row to clipboard:\n" + sage_code);		
	}
	else
	{	
	// Copy the text inside the text field
	var txt = table.cell( this ).data();
	navigator.clipboard.writeText(txt);

	// Alert the copied text
	alert("Copied to clipboard the text: " + txt + "\n (Use Ctrl-DblClick to copy Sage code for the row)");

    }});
    var table = $("#dgltable").DataTable({
	dom: 'lBfrtip',
	paging: false,
	ordering: false,
	fixedHeader: true,
	autoWidth: true,
	columnDefs: [ {
	    targets: [10, 11, 12],
	    render: $.fn.dataTable.render.ellipsis( 20 )
	}
        ],
	buttons: [
	    "searchBuilder",
	    "colvis",
	    "copyHtml5",
	    "csvHtml5"
	],
	initComplete: function () {
	    this.api()
		.columns()
		.every(function () {
		    let column = this;
		    let title = column.header().textContent;

		    if (title === "p" ||
			title === "type" ||
			title == 'n' ||
		        title === "label" ||
		        title === "char" ||
		        title === "trivial" ||
		        title === "recognized" ||
			title === "o(ζ)"
		       ) {

			// Create select element
			let select = document.createElement('select');
			select.add(new Option(''));
			column.header().append(select);

			// Apply listener for user change in value
			select.addEventListener('change', function () {
			    column
				.search(select.value, {exact: true})
				.draw();
			});

			// Add list of options
			column
			    .data()
			    .unique()
			    .sort()
			    .each(function (d, j) {
				select.add(new Option(d));
			    });

		    }
		    else
		    {
			// Create input element
			let input = document.createElement('input');
			input.placeholder = "Search " + title ;
			column.header().append(input);
			// Event listener for user input
			input.addEventListener('keyup', () => {
			    if (column.search() !== this.value) {
				column.search(input.value, true,false).draw();
			    }
			});

		    }
		});
	}
    });
	table.on('click', 'button', function (e) {
    let data = table.row(e.target.closest('tr')).data();
	let sage_code = "";
	let prec = Math.ceil((Math.log(data[12]) / Math.log(data[0])));
	if (prec == 0) {
		prec = 100;
	}
	if (data[5] == 'triv') {
		sage_code = "prime_list = [" + data[9].replace(/<[^>]*>?/gm, '').split("[")[1].split("]")[0] + "]; p = " + data[0] + "; prec = " + prec + "; label = " + data[1] + "; D = " + data[2] + "; n = " + data[3] + "; type = '" + data[4] + "'; char = '" + data[5] + "'; J0 = " + data[12] + "; J = Qq(p**2, prec, names='t')(2*J0);";
	}
	else {
		sage_code = "prime_list = [" + data[9].replace(/<[^>]*>?/gm, '').split("[")[1].split("]")[0] + "]; p = " + data[0] + "; prec = " + prec + "; label = " + data[1] + "; D = " + data[2] + "; n = " + data[3] + "; type = '" + data[4] + "'; char = '" + data[5] + "'; J0 = " + data[12] + "; J = Qq(p**2, prec, names='t')(J0) + ((1-Qq(p**2, prec, names='t')(J0)**2)/(t - t.trace()/2).norm()).sqrt() * (t - t.trace()/2);";
	}
	
	navigator.clipboard.writeText(sage_code);

	// Alert the copied text
	alert("Copied row to clipboard:\n" + sage_code);
});
});
;