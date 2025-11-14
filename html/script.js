
$(document).ready(function () {
	//Only needed for the filename of export files.
	//Normally set in the title tag of your page.
    document.title = "DGL Table";
    // DataTable initialisation
    $('#dgltable').on( 'dblclick', 'tbody td', function (e) {
	if (e.ctrlKey) {
    let data = table.row(e.target.closest('tr')).data();
	let sage_code = "";
	if (data[9] == '?') {
		sage_code = "prime_list = [];";
	}
	else
	{
		sage_code = "prime_list = ZZ(" + data[9].replace(/<[^>]*>?/gm, '').replaceAll('-','').split('=')[1].split('(')[0] + ").support(); ";
	}
	sage_code += " p = " + data[0] + "; prec = max(10, RR(" + data[12] + ").log(p).ceil());";
	sage_code += " label = '" + data[1] + "'; D = " + data[2] + "; n = " + data[3] + "; type = '" + data[4] + "'; char = '" + data[5] + "';";
	sage_code += " J0 = " + data[12] + "; Cp.<t> = Qq(p**2, prec);";
	if (data[5] == 'triv') {
		sage_code += "J = Cp(2*J0);";
	}
	else {
		sage_code += "J = Cp(J0) + ((1-Cp(J0)**2)/(t - t.trace()/2).norm()).sqrt() * (t - t.trace()/2);";
	}
	if (data[8] != '?') {
		sage_code += " L.<z> = NumberField(" + data[8] + "); "
	}
	if (data[11] != '?') {
		sage_code += " H.<ξ> = NumberField(" + data[11] + "); ";
	}
	
	navigator.clipboard.writeText(sage_code);

	// Alert the copied text
	// alert("Copied row to clipboard:\n" + sage_code);		
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
	buttons: [
	    "searchBuilder",
	    "colvis",
	    "copyHtml5",
	    "csvHtml5"
	],
	columnDefs: [ {
	    targets: [10, 11, 12],
	    render: $.fn.dataTable.render.ellipsis( 20 )
	},
	// Hide J by default
	{
		target: 12,
		visible: false
	},
	// Make discriminant column much narrower
	{
		target: 2,
		width: "120px"
	} 
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
			title === "o(ζ)"||
			title === "hE"
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
			    .sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
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
			// Minimize column width
			if (title === "D" || title === "J" || title === "hE" || title === "o(ζ)") {
			    input.style.width = "60px";
			}

		    }
		});
	}
    });
});
;

