
$(document).ready(function () {
	//Only needed for the filename of export files.
	//Normally set in the title tag of your page.
    document.title = "DGL Table";
    // DataTable initialisation
    $('#dgltable').on( 'dblclick', 'tbody td', function () {
	// Copy the text inside the text field
	var txt = table.cell( this ).data();
	navigator.clipboard.writeText(txt);

	// Alert the copied text
	alert("Copied to clipboard the text: " + txt);

    } );
    var table = $("#dgltable").DataTable({
	dom: 'lBfrtip',
	paging: false,
	ordering: false,
	fixedHeader: true,
	autoWidth: true,
	columnDefs: [ {
	    targets: [9, 10, 12],
	    render: $.fn.dataTable.render.ellipsis( 20 )
	} ],
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
});
