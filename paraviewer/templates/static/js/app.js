// ParaViewer Application Main JavaScript File

// Global variables
let currentViewer = null;
let variant_table;
let ndx;
let chromDimension;
let searchDimension;
let regionDimension;
let copyNumberDimension;
let familyIdDimension;
let paternalIdDimension;
let maternalIdDimension;
let sexDimension;
let phenotypeDimension;
let searchInput;
let regionInput;
let copyNumberInput;
let familyIdInput;
let paternalIdInput;
let maternalIdInput;
let sexInput;
let phenotypeInput;
let sizeChart;
let chromChart;
let variantCount;
let specialInfoFilterActive = false;
let showHelpActive = false;
let currentSexFilterValue = '';

// plot constraints
const plotw = 585;
const ploth = 150;

// Wait for DOM to be fully loaded before initializing
document.addEventListener('DOMContentLoaded', function () {
    // Configure DC.js colors
    dc.config.defaultColors(d3.schemeSet1);

    // Set up event handlers
    $('#filter-modal').on('hidden.bs.modal', function () {
        update_table();
    });

    // Initialize with the data variable from HTML template
    // The data variable is created by the jinja template in the HTML
    if (typeof window.appData !== 'undefined') {
        initializeApplication(window.appData);
    } else {
        console.error("Data not available");
    }
});

// Initialize the application
function initializeApplication(data) {
    try {
        ndx = crossfilter(data.data);  // Use data.data instead of just data
        var all = ndx.groupAll();

        // Initial dimension creation using the helper function
        const dims = createDimensions(ndx);
        chromDimension = dims.chrom;
        searchDimension = dims.search;
        regionDimension = dims.region;
        copyNumberDimension = dims.copyNumber;
        familyIdDimension = dims.familyId;
        paternalIdDimension = dims.paternalId;
        maternalIdDimension = dims.maternalId;
        sexDimension = dims.sex;
        phenotypeDimension = dims.phenotype;

        // Initialize table filters
        searchInput = dc.textFilterWidget("#sample-search");
        regionInput = dc.textFilterWidget("#region-search");
        copyNumberInput = dc.textFilterWidget("#copy-number-search");
        
        // Only initialize pedigree column filters if they exist in the DOM
        if (document.getElementById("family-id-search")) {
            familyIdInput = dc.textFilterWidget("#family-id-search");
        }
        if (document.getElementById("paternal-id-search")) {
            paternalIdInput = dc.textFilterWidget("#paternal-id-search");
        }
        if (document.getElementById("maternal-id-search")) {
            maternalIdInput = dc.textFilterWidget("#maternal-id-search");
        }
        // Sex uses a standalone select, not a DC widget
        if (document.getElementById("phenotype-search")) {
            phenotypeInput = dc.textFilterWidget("#phenotype-search");
        }
        chromChart = dc.barChart("#chrom-chart");

        build_table(chromDimension.top(Infinity));
        var chromGroup = chromDimension.group().reduceCount();
        var nonEmptyChromGroup = remove_empty_bins(chromGroup);

        setupCharts(ndx, all, nonEmptyChromGroup);

        dc.renderAll();
        
        // Populate Sex dropdown when the filter modal opens
        $('#filter-modal').off('shown.bs.modal.setupSex').on('shown.bs.modal.setupSex', function () {
            if (document.getElementById("sex-search")) {
                setupSexDropdown();
            }
        });

    } catch (error) {
        console.error("Error initializing application:", error);
    }
}

// Set up charts
function setupCharts(ndx, all, nonEmptyChromGroup) {
    // Helper function to set up text filter widgets
    function setupTextFilter(input, dimension) {
        input
            .dimension(dimension)
            .on('renderlet', function () {
                d3.selectAll(".dc-text-filter-input")
                    .classed("form-control", true);
                d3.selectAll(`#${input.anchorName()}.dc-chart`)
                    .classed("col-12", true);
            });
    }

    // Set up all text filters
    setupTextFilter(searchInput, searchDimension);
    setupTextFilter(regionInput, regionDimension);
    setupTextFilter(copyNumberInput, copyNumberDimension);
    
    // Only set up pedigree column filters if they were initialized
    if (familyIdInput) {
        setupTextFilter(familyIdInput, familyIdDimension);
    }
    if (paternalIdInput) {
        setupTextFilter(paternalIdInput, paternalIdDimension);
    }
    if (maternalIdInput) {
        setupTextFilter(maternalIdInput, maternalIdDimension);
    }
    // Sex handled via standalone select; no DC widget setup here
    if (phenotypeInput) {
        setupTextFilter(phenotypeInput, phenotypeDimension);
    }

    // chromosome
    chromChart
        .width(plotw).height(ploth).gap(1)
        .margins({ top: 10, right: 50, bottom: 70, left: 40 })
        .x(d3.scaleBand())
        .xUnits(dc.units.ordinal)
        .yAxisLabel('Count')
        .elasticX(true)
        .elasticY(true)
        .dimension(chromDimension)
        .group(nonEmptyChromGroup)
        .ordering((d) => {
            // Strip out 'chr' prefix if it exists
            let key = d.key.toLowerCase().replace('chr', '');
            let v = parseInt(key);
            if (v) {
                return v;
            } else {
                // Handle special cases like 'X', 'Y', 'M'
                if (key === 'x') return 23;
                if (key === 'y') return 24;
                if (key === 'm' || key === 'mt') return 25;
                return key;
            }
        })
        .on('renderlet', function (chart) {
            chart.selectAll('g.x text')
                .attr('transform', 'rotate(90)')
                .attr('text-anchor', 'start')
                .attr('y', -5)
                .attr('x', 10);
        });
    chromChart.yAxis().ticks(5);
}

function setupSexDropdown() {
    // Build unique values
    const allData = ndx.all();
    const uniqueSexValues = [...new Set(allData.map(d => d.Sex).filter(sex => sex && sex !== ''))].sort();

    const select = document.getElementById('sex-search');
    if (!select) return;

    // Determine current filter value
    // Prefer the explicitly tracked value; fall back to deriving from dimension once
    let derivedValue = '';
    try {
        if (sexDimension && typeof sexDimension.hasFilter === 'function' && sexDimension.hasFilter()) {
            const filteredData = sexDimension.top(2);
            if (filteredData && filteredData.length > 0) {
                const vals = [...new Set(filteredData.map(d => d.Sex).filter(Boolean))];
                if (vals.length === 1) derivedValue = vals[0];
            }
        }
    } catch (_) { /* ignore */ }

    // Populate
    select.innerHTML = '';
    const placeholder = document.createElement('option');
    placeholder.value = '';
    placeholder.textContent = 'Select Sex...';
    select.appendChild(placeholder);

    uniqueSexValues.forEach(sex => {
        const option = document.createElement('option');
        option.value = sex;
        option.textContent = sex;
        select.appendChild(option);
    });

    // Set the current filter value if one exists
    const valueToSet = currentSexFilterValue || derivedValue || '';
    if (valueToSet) {
        if (uniqueSexValues.includes(valueToSet)) {
            select.value = valueToSet;
        } else {
            // Case-insensitive fallback
            const match = uniqueSexValues.find(v => String(v).toLowerCase() === String(valueToSet).toLowerCase());
            if (match) {
                select.value = match;
            } else {
                select.value = '';
            }
        }
    } else {
        select.value = '';
    }

    // Wire filtering (overwrite to avoid duplicate listeners)
    select.onchange = function () {
        const val = this.value;
        if (!val) {
            sexDimension.filter(null);
        } else {
            sexDimension.filter(val);
        }
        currentSexFilterValue = val || '';
        dc.redrawAll();
        update_table();
    };
}

function handleKeyPress(event) {
    // If no viewer is active, don't do anything
    if (!currentViewer) return;

    // Prevent default scrolling behavior for left/right keys
    if (event.key === 'ArrowLeft' || event.key === 'ArrowRight') {
        event.preventDefault();

        const current = $('tr.selected');
        if (!current.length) return;  // No row is selected

        const prev = current.prev();
        const next = current.next();

        try {
            if (event.key === 'ArrowLeft' && prev.length > 0) {
                currentViewer.destroy();
                currentViewer = null;
                table_click(prev[0], variant_table);
            } else if (event.key === 'ArrowRight' && next.length > 0) {
                currentViewer.destroy();
                currentViewer = null;
                table_click(next[0], variant_table);
            }
        } catch (error) {
            console.error("Error handling key navigation:", error);
            if (currentViewer) {
                try {
                    currentViewer.destroy();
                } catch (e) { /* ignore error during cleanup */ }
                currentViewer = null;
            }
            document.body.classList.remove('viewer-open');
        }
    }
}

function showImageFromHash() {
    if (window.location.hash) {
        try {
            // Format is #chrom:start-end:sample
            const hash = window.location.hash.substring(1);
            const match = hash.match(/^([^:]+):(\d+)-(\d+):(.+)$/);

            if (match) {
                const [_, chrom, startStr, endStr, sample] = match;
                const start = parseInt(startStr);
                const end = parseInt(endStr);

                if (isNaN(start) || isNaN(end)) {
                    console.warn(`Invalid coordinates in hash: ${hash}`);
                    history.pushState(null, '', window.location.pathname);
                    return;
                }

                // Find the matching row
                let foundRow = false;
                let matchingRowIndex = -1;

                variant_table.rows().every(function (rowIndex) {
                    const data = this.data();

                    if (data.Chrom === chrom &&
                        data.Start === start &&
                        data.End === end &&
                        data.Sample === sample) {
                        matchingRowIndex = rowIndex;
                        foundRow = true;
                        return false; // break the loop
                    }
                    return true;
                });

                // Handle case where a matching row was found - use the index directly
                if (foundRow && matchingRowIndex >= 0) {
                    // Use the row index directly - no need to find the DOM node
                    table_click(matchingRowIndex, variant_table);
                } else {
                    console.warn(`No matching row found for ${chrom}:${start}-${end}:${sample}`);
                    // Clear the hash to prevent repeated errors
                    history.pushState(null, '', window.location.pathname);
                }
            } else {
                console.warn(`Hash does not match expected format: ${hash}`);
                history.pushState(null, '', window.location.pathname);
            }
        } catch (error) {
            console.error("Error processing hash:", error);
            history.pushState(null, '', window.location.pathname);
        }
    }
}

// Function to handle Quick view button clicks
function quickViewImage(chrom, start, end, sample, imagePath, event) {
    event.preventDefault();
    event.stopPropagation();
    
    // Find the row index for this data
    let rowIndex = -1;
    variant_table.rows().every(function (index) {
        const rowData = this.data();
        if (rowData.Chrom === chrom && 
            rowData.Start === start && 
            rowData.End === end && 
            rowData.Sample === sample) {
            rowIndex = index;
            return false; // break the loop
        }
        return true;
    });
    
    if (rowIndex >= 0) {
        table_click(rowIndex, variant_table);
    } else {
        console.warn("Could not find matching row for Quick view");
    }
}

const table_click = (selection, table) => {
    // Ensure we have a valid selection before proceeding
    if (!selection || !table) {
        console.warn("Invalid selection");
        return;
    }

    let rowIndex;

    // Handle both DOM nodes and direct row indices
    if (typeof selection === 'number') {
        // It's already a row index
        rowIndex = selection;
    } else {
        // It's a DOM node
        rowIndex = table.row(selection).index();
    }

    if (rowIndex === undefined) {
        console.warn("Could not determine row index");
        return;
    }

    // Get row data before making any DOM changes
    const row = table.row(rowIndex).data();
    if (!row) {
        console.warn("Selected row has no data");
        return;
    }

    // Check if we have a valid row with the required properties
    if (!row.Chrom || !row.Start || !row.End || !row.Sample || !row.Image) {
        console.warn("Selected row is missing required data");
        return;
    }

    // Try to get the DOM node for visual selection
    const rowNode = table.row(rowIndex).node();

    // Only update visual selection if we have a valid DOM node
    if (rowNode) {
        table.$('tr.selected').removeClass('selected');
        $(rowNode).addClass('selected');
    }

    // Get references to adjacent rows for navigation buttons
    const currentRowIndex = rowIndex;
    const hasPrevRow = currentRowIndex > 0;
    const hasNextRow = currentRowIndex < table.rows().count() - 1;

    // Update URL hash with row information including sample
    window.history.pushState(null, '',
        `#${row.Chrom}:${row.Start}-${row.End}:${row.Sample}`);

    // Create a new Image object
    let img = new Image();
    img.src = row.Image;

    img.onerror = function () {
        console.error(`Image not found: ${img.src}`);
    };

    img.onload = function () {
        // Image loaded successfully
    };

    // Create and configure the viewer
    let viewer = new Viewer(img, {
        hidden: function () {
            document.removeEventListener('keydown', handleKeyPress);
            currentViewer = null;
            viewer.destroy();

            // Re-enable default browser keyboard behavior when viewer is closed
            document.body.classList.remove('viewer-open');
        },
        title: function () {
            return `${row.Sample} - ${row.Region} - ${row.Chrom}:${row.Start}-${row.End}`;
        },
        shown: function () {
            // Disable default browser keyboard behavior when viewer is open
            document.body.classList.add('viewer-open');
        },
        toolbar: {
            zoomIn: 4,
            zoomOut: 4,
            oneToOne: 4,
            reset: 4,
            prev: {
                show: hasNextRow,  // Show "prev" button when there's a next row
                size: "large",
                click: function () {
                    try {
                        if (!hasNextRow) return;

                        viewer.destroy();
                        currentViewer = null;

                        // Get the next row using DataTables API and DOM element
                        const nextRow = $(variant_table.row(currentRowIndex + 1).node());
                        if (nextRow.length) {
                            table_click(nextRow[0], variant_table);
                        }
                    } catch (error) {
                        console.error("Error navigating to next image:", error);
                        // Try to clean up
                        try {
                            viewer.destroy();
                        } catch (e) { /* ignore error during cleanup */ }
                        currentViewer = null;
                    }
                }
            },
            play: { show: false },
            next: {
                show: hasPrevRow,  // Show "next" button when there's a previous row
                size: "large",
                click: function () {
                    try {
                        if (!hasPrevRow) return;

                        viewer.destroy();
                        currentViewer = null;

                        // Get the previous row using DataTables API and DOM element
                        const prevRow = $(variant_table.row(currentRowIndex - 1).node());
                        if (prevRow.length) {
                            table_click(prevRow[0], variant_table);
                        }
                    } catch (error) {
                        console.error("Error navigating to previous image:", error);
                        // Try to clean up
                        try {
                            viewer.destroy();
                        } catch (e) { /* ignore error during cleanup */ }
                        currentViewer = null;
                    }
                }
            },
            rotateLeft: { show: false },
            rotateRight: { show: false },
            flipHorizontal: { show: false },
            flipVertical: { show: false },
        },
        transition: false,
        navbar: false,
    });

    // Add the event listener when viewer is shown
    currentViewer = viewer;
    document.addEventListener('keydown', handleKeyPress);

    // Show the viewer
    viewer.show();
}

// Function to create filterable column definitions
function createFilterableColumn(index, filterFn = 'filterByColumn') {
    return {
        targets: index,
        render: function (data, type, row) {
            if (type === 'display' && data != null) {
                return `<a href="#" onclick="return filterByColumn('${data}', ${index}, event)" style="color: inherit; text-decoration: none;">${data}</a>`;
            }
            return data;
        }
    };
}

function build_table(data) {
    // hide the placeholder and show the datatable
    d3.select('#variant-table-placeholder').property("hidden", true)
    d3.select('#variant-table-div').property("hidden", false)

    // Check if there are any trio records
    const hasTrioRecords = data.some(row => 
        row.PaternalID && row.PaternalID !== 'NA' && 
        row.MaternalID && row.MaternalID !== 'NA'
    );

    let cols = [
        { data: null, title: 'Quick view', defaultContent: '' },
        { data: 'Chrom', title: 'Chrom' },
        { data: 'Start', title: 'Start' },
        { data: 'End', title: 'End' },
        { data: 'Region', title: 'Region' },
        { data: 'Sample', title: 'Sample' },
        { data: 'CopyNumber', title: 'Copy Number' },
        { data: 'FamilyID', title: 'Family ID', visible: window.appData.has_pedigree_columns.FamilyID },
        { data: 'PaternalID', title: 'Paternal ID', visible: window.appData.has_pedigree_columns.PaternalID },
        { data: 'MaternalID', title: 'Maternal ID', visible: window.appData.has_pedigree_columns.MaternalID },
        { data: 'Sex', title: 'Sex', visible: window.appData.has_pedigree_columns.Sex },
        { data: 'Phenotype', title: 'Phenotype', visible: window.appData.has_pedigree_columns.Phenotype },
        { data: null, title: 'View Trio', defaultContent: '', visible: hasTrioRecords },
        { data: null, title: 'Open IGV', defaultContent: '' },
        { data: null, title: 'IGV Session', defaultContent: '' }
    ]

    variant_table = $("#variant-table").DataTable({
        data: data,
        columns: cols,
        deferRender: true,
        scrollY: '80vh',
        scrollCollapse: true,
        scroller: {
            loadingIndicator: true,
            displayBuffer: 20,
            serverWait: 100
        },
        info: true,
        buttons: [
            {
                extend: 'copyHtml5',
                fieldSeparator: '\t',
                fieldBoundary: '',
                exportOptions: {
                    columns: [1, 2, 3, 4, 5, 6],
                    format: {
                        body: function (data, row, column, node) {
                            const rowData = variant_table.row(row).data();
                            if (column === 5) {
                                return rowData.sample;
                            }
                            return data;
                        }
                    }
                }
            },
            {
                extend: 'csvHtml5',
                fieldSeparator: '\t',
                fieldBoundary: '',
                extension: '.tsv',
                filename: 'paraphase_variants',
                exportOptions: {
                    columns: [1, 2, 3, 4, 5, 6],
                    format: {
                        body: function (data, row, column, node) {
                            const rowData = variant_table.row(row).data();
                            if (column === 5) {
                                return rowData.sample;
                            }
                            return data;
                        }
                    }
                }
            }
        ],
        dom: 'Brti',
        initComplete: function () {
            $('.dt-buttons').hide();
            
            // Synchronize header and body scroll for horizontal scrolling
            const scrollBody = document.querySelector('#variant-table_wrapper .dataTables_scrollBody');
            const scrollHead = document.querySelector('#variant-table_wrapper .dataTables_scrollHead');
            
            if (scrollBody && scrollHead) {
                // When body scrolls, update header
                scrollBody.addEventListener('scroll', function() {
                    scrollHead.scrollLeft = scrollBody.scrollLeft;
                });
                
                // When header scrolls, update body
                scrollHead.addEventListener('scroll', function() {
                    scrollBody.scrollLeft = scrollHead.scrollLeft;
                });
            }
            
            setTimeout(function () {
                showImageFromHash();
            }, 200);
        },
        infoCallback: (oSettings, iStart, iEnd, iMax, iTotal, sPre) => {
            return `
        <span class="datatable-info"> 
            <span class="pr-2">Showing <b>${iStart}</b> - <b>${iEnd}</b> of <b>${iTotal}</b> records</span>
            <button type="button" class="btn btn-outline-primary btn-sm" data-toggle="modal" data-target="#filter-modal" title="Show filters">
                <span class="fas fa-filter"></span>
            </button>
            <button type="button" class="btn btn-outline-primary btn-sm mr-2" onclick="resetFilter()" title="Reset all filters">
                <span class="fas fa-undo"></span>
            </button>
            <span class="dropup">
                <button type="button" class="btn btn-sm btn-outline-primary dropdown-toggle" id="download-menu" title="Save table" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                    <span class="fas fa-save"></span>
                </button>
                <span class="dropdown-menu" aria-labelledby="download-menu">
                    <h6 class="dropdown-header">Save ${iTotal} rows as:</h6>
                    <button class="dropdown-item" type="button" id="tsv-button-download" onclick="tsv_button_click()">
                        TSV
                    </button>
                    <button class="dropdown-item" type="button" id="copy-button-download" onclick="copy_button_click()">
                        Copy
                    </button>
                </span>
            </span>
            <div class="form-check form-check-inline ml-2">
                <input class="form-check-input" type="checkbox" id="special-info-filter" onchange="toggleSpecialInfoFilter(this)" ${specialInfoFilterActive ? 'checked' : ''}>
                <label class="form-check-label" for="special-info-filter">Special Info Only</label>
            </div>
            <div class="form-check form-check-inline ml-2">
                <input class="form-check-input" type="checkbox" id="show-help" onchange="toggleShowHelp(this)" ${showHelpActive ? 'checked' : ''}>
                <label class="form-check-label" for="show-help">Show Help</label>
            </div>
        </span>
        `
        },
        columnDefs: [
            // Quick view button column
            {
                targets: 0,
                render: function (data, type, row) {
                    if (type === 'display') {
                        const infoIndicator = row.SpecialInfo ?
                            '<span class="info-indicator" title="Has additional information">★</span>' :
                            '<span class="info-placeholder"></span>';
                        return `<div class="quick-view-container">
                                    ${infoIndicator}
                                    <button class="btn btn-sm btn-outline-primary quick-view-btn" 
                                        onclick="quickViewImage('${row.Chrom}', ${row.Start}, ${row.End}, '${row.Sample}', '${row.Image}', event)"
                                        title="Quick view image">
                                        <i class="fas fa-eye"></i>
                                    </button>
                                </div>`;
                    }
                    return '';
                },
                className: 'text-center',
                orderable: false,
                searchable: false
            },
            // Special case for Chrom column
            {
                targets: 1,
                render: function (data, type, row) {
                    if (type === 'display' && data != null) {
                        return `<a href="#" onclick="filterByColumn('${data}', 1, event)" style="color: inherit; text-decoration: none;">${data}</a>`;
                    }
                    return data;
                }
            },
            // Common columns using helper
            createFilterableColumn(2),  // Start
            createFilterableColumn(3),  // End
            createFilterableColumn(4, 'filterByRegion'),  // Region
            createFilterableColumn(5, 'filterBySample'),  // Sample
            createFilterableColumn(6),  // Copy Number
            createFilterableColumn(7),  // Family ID
            createFilterableColumn(8),  // Paternal ID
            createFilterableColumn(9),  // Maternal ID
            createFilterableColumn(10), // Sex
            createFilterableColumn(11), // Phenotype
            // Special case for View Trio column
            {
                targets: 12,
                render: function (data, type, row) {
                    if (type === 'display') {
                        // Get all unique sample IDs from the data
                        const allSamples = new Set(window.appData.data.map(d => d.Sample));

                        // Only show button if both parents are known and exist as samples
                        if (row.PaternalID && row.PaternalID !== 'NA' &&
                            row.MaternalID && row.MaternalID !== 'NA' &&
                            allSamples.has(row.PaternalID) && allSamples.has(row.MaternalID)) {
                            return `<button class="btn btn-sm btn-outline-primary view-trio-btn" 
                                        onclick="viewTrio('${row.Sample}', '${row.PaternalID}', '${row.MaternalID}', event)"
                                        title="View sample and parents">
                                        <i class="fas fa-users"></i>
                                    </button>`;
                        }
                        return '';
                    }
                    return '';
                },
                className: 'text-center',
                orderable: false,
                searchable: false
            },
            // Special case for IGV columns
            {
                targets: 13,  // Open IGV column
                render: function (data, type, row) {
                    if (type === 'display') {
                        const isLocalFile = window.location.protocol === 'file:';
                        const titleText = isLocalFile ?
                            'Select project directory to view IGV (required due to browser security restrictions)' :
                            'Open IGV with BAM track';

                        return `<button class="btn btn-sm btn-outline-primary" 
                                    onclick="openIGV('${row.Chrom}', ${row.Start}, ${row.End}, '${row.BAM}', event)"
                                    title="${titleText}">
                                    <i class="fas fa-external-link-alt"></i>
                                </button>`;
                    }
                    return '';
                },
                className: 'text-center',
                orderable: false,
                searchable: false
            },
            {
                targets: 14,  // IGV Session column
                render: function (data, type, row) {
                    if (type === 'display') {
                        if (row.IGVSession) {
                            const isTrio = row.Sample.endsWith('-trio');
                            if (isTrio) {
                                return `<button class="btn btn-sm btn-outline-primary igv-download-btn" 
                                            onclick="downloadFile('${row.IGVSession}', event)" 
                                            title="Download IGV Session File">
                                            <i class="fas fa-download"></i> IGV
                                        </button>`;
                            } else {
                                return `<div class="btn-group btn-group-sm">
                                            <button class="btn btn-sm btn-outline-primary igv-download-btn" 
                                                onclick="downloadFile('${row.IGVSession}', event)" 
                                                title="Download IGV Session File">
                                                <i class="fas fa-download"></i> IGV
                                            </button>
                                            <button class="btn btn-sm btn-outline-primary" 
                                                onclick="downloadFile('${row.BAM}', event)" 
                                                title="Download BAM File">
                                                <i class="fas fa-download"></i> BAM
                                            </button>
                                            <button class="btn btn-sm btn-outline-primary" 
                                                onclick="downloadFile('${row.BAI}', event)" 
                                                title="Download BAI File">
                                                <i class="fas fa-download"></i> BAI
                                            </button>
                                        </div>`;
                            }
                        }
                        return '';
                    }
                    return '';
                },
                className: 'text-center',
                orderable: false,
                searchable: false
            }
        ],
        lengthChange: false,
        order: [[1, 'asc'], [2, 'asc']],
        rowId: function (data) {
            return `row-${data.Chrom}-${data.Start}-${data.End}-${data.Sample.replace(/\s+/g, '-')}`;
        },
        createdRow: function (row, data, dataIndex) {
            if (data.SpecialInfo) {
                $(row).addClass('has-info-dropdown');

                $(row).on('mousemove', function (e) {
                    $('#info-popup').remove();

                    const dropdownContent = document.createElement('div');
                    dropdownContent.id = 'info-popup';
                    dropdownContent.className = 'info-dropdown-content';
                    dropdownContent.style.maxWidth = '600px';  // Double the default max-width
                    dropdownContent.style.whiteSpace = 'normal';  // Allow text to wrap
                    dropdownContent.style.wordBreak = 'break-word';  // Break long words if needed
                    dropdownContent.innerHTML = `
                        <h6>Special Info:</h6>
                        <table class="special-info-table">
                            ${data.SpecialInfo.split(';').filter(item => item.trim()).map(row => {
                        const cells = row.split(',').map(cell => cell.trim());
                        return `<tr>${cells.map(cell => `<td>${cell}</td>`).join('')}</tr>`;
                    }).join('')}
                        </table>
                    `;

                    document.body.appendChild(dropdownContent);

                    const dropdown = $('#info-popup');
                    
                    // Smart positioning function
                    const positionModal = (element, x, y) => {
                        element.css({ display: 'block', left: '-9999px', top: '-9999px' });
                        const rect = element[0].getBoundingClientRect();
                        const viewport = { w: window.innerWidth, h: window.innerHeight };
                        
                        const left = x + rect.width > viewport.w ? x - rect.width - 10 : x + 10;
                        const top = y + rect.height > viewport.h ? y - rect.height - 10 : y + 10;
                        
                        element.css({ left: Math.max(10, left) + 'px', top: Math.max(10, top) + 'px' });
                    };
                    
                    positionModal(dropdown, e.pageX, e.pageY);
                });

                $(row).on('mouseleave', function () {
                    $('#info-popup').remove();
                });
            }
        }
    });

    // Row click handler removed - now using Quick view buttons instead

    // Handle context menu (right-click)
    setupContextMenu();
}

function tsv_button_click() {
    // Get all data that passes the current filters
    let filteredData = [];
    variant_table.rows({ search: 'applied' }).every(function () {
        const rowData = this.data();
        filteredData.push([
            rowData.Chrom,
            rowData.Start,
            rowData.End,
            rowData.Region,
            rowData.Sample,
            rowData.CopyNumber,
            rowData.FamilyID || '',
            rowData.PaternalID || '',
            rowData.MaternalID || '',
            rowData.Sex || '',
            rowData.Phenotype || '',
            rowData.SpecialInfo || ''
        ]);
    });

    // Create TSV content (excluding Quick view column)
    let tsvContent = "Chrom\tStart\tEnd\tRegion\tSample\tCopy Number\tFamily ID\tPaternal ID\tMaternal ID\tSex\tPhenotype\tSpecial Info\n";
    filteredData.forEach(row => {
        tsvContent += row.join('\t') + '\n';
    });

    // Create and trigger download
    const blob = new Blob([tsvContent], { type: 'text/tab-separated-values' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.style.display = 'none';
    a.href = url;
    a.download = 'paraphase_variants.tsv';
    document.body.appendChild(a);
    a.click();
    window.URL.revokeObjectURL(url);
    document.body.removeChild(a);
}

function copy_button_click() {
    // Get all data that passes the current filters
    let filteredData = [];
    variant_table.rows({ search: 'applied' }).every(function () {
        const rowData = this.data();
        filteredData.push([
            rowData.Chrom,
            rowData.Start,
            rowData.End,
            rowData.Region,
            rowData.Sample,
            rowData.CopyNumber,
            rowData.FamilyID || '',
            rowData.PaternalID || '',
            rowData.MaternalID || '',
            rowData.Sex || '',
            rowData.Phenotype || '',
            rowData.SpecialInfo || ''
        ]);
    });

    // Create TSV content for clipboard (excluding Quick view column)
    let tsvContent = "Chrom\tStart\tEnd\tRegion\tSample\tCopy Number\tFamily ID\tPaternal ID\tMaternal ID\tSex\tPhenotype\tSpecial Info\n";
    filteredData.forEach(row => {
        tsvContent += row.join('\t') + '\n';
    });

    // Copy to clipboard
    const textarea = document.createElement('textarea');
    textarea.value = tsvContent;
    document.body.appendChild(textarea);
    textarea.select();
    document.execCommand('copy');
    document.body.removeChild(textarea);

    // Show a notification that data was copied
    alert(`Copied ${filteredData.length} rows to clipboard`);
}

function update_table() {
    variant_table.clear().rows.add(chromDimension.top(Infinity)).draw();
}

function remove_empty_bins(source_group) {
    return {
        all: function () {
            return source_group.all().filter(function (d) {
                return d.value != 0
            })
        }
    }
}

// https://jsfiddle.net/gordonwoodhull/g34Ldwaz/8/
// https://github.com/dc-js/dc.js/issues/348
function index_group(group) {
    return {
        all: function () {
            return group.all().map(function (kv, i) {
                return { key: i, value: kv.value }
            })
        }
    }
}

// Function to create all dimensions
function createDimensions(ndx) {
    return {
        chrom: ndx.dimension(d => d.Chrom),
        search: ndx.dimension(d => d.Sample),
        region: ndx.dimension(d => d.Region),
        copyNumber: ndx.dimension(d => String(d.CopyNumber)),
        familyId: ndx.dimension(d => d.FamilyID),
        paternalId: ndx.dimension(d => d.PaternalID),
        maternalId: ndx.dimension(d => d.MaternalID),
        sex: ndx.dimension(d => d.Sex),
        phenotype: ndx.dimension(d => String(d.Phenotype))
    };
}

// Helper function to map column indices to filter input selectors
function getFilterInputSelector(columnIndex) {
    // Check if pedigree columns are available
    const hasPedigreeColumns = window.appData && window.appData.has_pedigree_columns;
    
    const columnToInputMap = {
        1: "#sample-search",      // Chrom -> Sample (no specific filter for Chrom)
        2: "#sample-search",      // Start -> Sample (no specific filter for Start)
        3: "#sample-search",      // End -> Sample (no specific filter for End)
        4: "#region-search",      // Region
        5: "#sample-search",      // Sample
        6: "#copy-number-search", // Copy Number
        7: hasPedigreeColumns && hasPedigreeColumns.FamilyID ? "#family-id-search" : null,   // Family ID
        8: hasPedigreeColumns && hasPedigreeColumns.PaternalID ? "#paternal-id-search" : null, // Paternal ID
        9: hasPedigreeColumns && hasPedigreeColumns.MaternalID ? "#maternal-id-search" : null, // Maternal ID
        10: hasPedigreeColumns && hasPedigreeColumns.Sex ? "#sex-search" : null,        // Sex
        11: hasPedigreeColumns && hasPedigreeColumns.Phenotype ? "#phenotype-search" : null  // Phenotype
    };
    return columnToInputMap[columnIndex];
}

// Consolidated filter functions
function filterByColumn(value, columnIndex, event) {
    event.preventDefault();
    event.stopPropagation();

    // Handle Sex via dropdown + Crossfilter exact filter
    const inputSelector = getFilterInputSelector(columnIndex);
    if (inputSelector === '#sex-search') {
        const select = document.getElementById('sex-search');
        if (select) select.value = value;
        if (!value) {
            sexDimension.filter(null);
        } else {
            sexDimension.filter(value);
        }
        currentSexFilterValue = value || '';
        dc.redrawAll();
        update_table();
        return false;
    }

    // Default behavior for other columns: DataTables exact match of clicked value
    const searchTerms = value.toLowerCase().split(/\s+/).filter(term => term.length > 0);

    $.fn.dataTable.ext.search.push((settings, data) => {
        const sampleValue = (data[columnIndex] || '').toLowerCase();
        return searchTerms.some(term => sampleValue === term);
    });

    // Update the corresponding input field in the filter modal for text widgets
    if (inputSelector && inputSelector !== '#sex-search') {
        const input = document.querySelector(`${inputSelector} .dc-text-filter-input`);
        if (input) input.value = value;
    }

    variant_table.draw();
}

// Setup context menu functionality
function setupContextMenu() {
    const menu = document.getElementById('context-menu');
    const table = document.querySelector('#variant-table');
    const copyPageUrlItem = document.getElementById('copy-page-url');
    const viewImageItem = document.getElementById('view-image');

    table.addEventListener('contextmenu', e => {
        e.preventDefault();
        const tr = $(e.target).closest('tr')[0];
        if (!tr) return;

        menu.style.cssText = `top: ${e.pageY}px; left: ${e.pageX}px; display: block;`;
        menu.dataset.rowId = tr.id;
    });

    document.addEventListener('click', () => menu.style.display = 'none');

    // Copy Page URL with hash
    copyPageUrlItem.addEventListener('click', () => {
        const tr = document.getElementById(menu.dataset.rowId);
        if (!tr) return;

        const rowData = variant_table.row(tr).data();
        if (!rowData || !rowData.Chrom || !rowData.Start || !rowData.End || !rowData.Sample) return;

        const hash = `#${rowData.Chrom}:${rowData.Start}-${rowData.End}:${rowData.Sample}`;
        const pageUrl = `${window.location.origin}${window.location.pathname}${hash}`;
        copyToClipboard(pageUrl);
    });

    // View image menu item
    viewImageItem.addEventListener('click', () => {
        try {
            const tr = document.getElementById(menu.dataset.rowId);
            if (!tr) return;
            table_click(tr, variant_table);
        } catch (error) {
            console.error("Error handling view image action:", error);
        }
    });
}

function resetFilter() {
    // 1) Clear DC text filter widgets and their DOM inputs BEFORE redraw
    //    so widgets don't reapply filters during render
    const widgetVars = [
        searchInput,
        regionInput,
        copyNumberInput,
        familyIdInput,
        paternalIdInput,
        maternalIdInput,
        phenotypeInput
    ].filter(Boolean);
    widgetVars.forEach(widget => {
        try { widget.filter(null); } catch (e) { /* no-op */ }
    });

    const textSelectors = [
        "#sample-search",
        "#region-search",
        "#copy-number-search",
        "#family-id-search",
        "#paternal-id-search",
        "#maternal-id-search",
        "#phenotype-search"
    ];
    textSelectors.forEach(selector => {
        const el = document.querySelector(`${selector} .dc-text-filter-input`);
        if (el) {
            el.value = '';
            // Ensure widget listeners see the change
            try { el.dispatchEvent(new Event('input', { bubbles: true })); } catch (e) {}
        }
    });

    // 2) Reset Sex dropdown and clear filter
    const sexSelect = document.getElementById('sex-search');
    if (sexSelect) sexSelect.value = '';
    if (sexDimension) sexDimension.filter(null);
    currentSexFilterValue = '';

    // 3) Clear any DataTables search filters and custom search fns
    variant_table.search('').columns().search('');
    while ($.fn.dataTable.ext.search.length > 0) {
        $.fn.dataTable.ext.search.pop();
    }

    // Reset Special Info Only checkbox UI and state
    specialInfoFilterActive = false;
    const specialInfoCheckbox = document.getElementById('special-info-filter');
    if (specialInfoCheckbox) specialInfoCheckbox.checked = false;

    // 4) Clear Crossfilter/DC filters and redraw
    dc.filterAll();
    dc.redrawAll();

    // 5) Refresh the table from Crossfilter
    update_table();

    // 6) Clear the URL hash
    history.pushState(null, '', window.location.pathname);
}

/**
 * Download a file (IGV session, BAM, or BAI)
 * @param {string} filePath - Path to the file
 * @param {Event} event - The click event
 */
function downloadFile(filePath, event) {
    event.stopPropagation();
    const fileName = filePath.split('/').pop() || 'file';
    const isXmlFile = fileName.toLowerCase().endsWith('.xml');

    if (isXmlFile) {
        window.open(filePath, '_blank');
        return;
    }

    const a = document.createElement('a');
    a.href = filePath;
    a.download = fileName;
    a.rel = 'noopener noreferrer';
    document.body.appendChild(a).click();
    document.body.removeChild(a);
}

// Consolidated event handlers
function setupEventHandlers() {
    document.addEventListener('click', (e) => {
        if (e.target.matches('.reset-btn')) resetAll();
        if (e.target.matches('.export-btn')) handleExport(e.target.dataset.format);
    });
}

// Simplified URL parameter handling
function updateUrlParams() {
    const params = new URLSearchParams();
    if (searchInput.value()) params.set('search', searchInput.value());
    if (chromChart.filters().length) params.set('chrom', chromChart.filters().join(','));
    history.replaceState(null, '', `?${params.toString()}`);
}

// Streamlined search function
function performSearch(text) {
    searchDimension.filter(text);
    updateUrlParams();
    dc.redrawAll();
    update_table();
}

// Consolidated loading functions
function toggleLoading(show) {
    document.getElementById('loading-spinner').style.display = show ? 'block' : 'none';
}

// Simplified error handling
function showError(message) {
    const errorDiv = document.getElementById('error-message');
    errorDiv.textContent = message;
    errorDiv.style.display = 'block';
    setTimeout(() => errorDiv.style.display = 'none', 5000);
}

// Restore clipboard functionality
function copyToClipboard(text) {
    const textarea = document.createElement('textarea');
    textarea.value = text;
    document.body.appendChild(textarea);
    textarea.select();

    try {
        const successful = document.execCommand('copy');
        if (successful) {
            const successMsg = document.querySelector('.clipboard-success');
            successMsg.style.display = 'block';
            setTimeout(() => successMsg.style.display = 'none', 2000);
        }
    } catch (err) {
        console.error('Unable to copy to clipboard', err);
        alert('Failed to copy URL to clipboard');
    }

    document.body.removeChild(textarea);
}

// Add new function to handle IGV opening
function openIGV(chrom, start, end, bamPath, event) {
    event.stopPropagation();
    const isLocalFile = window.location.protocol === 'file:';

    // Get the current row data to access BAM and BAI paths
    const currentRow = variant_table.row($(event.target).closest('tr')).data();
    if (!currentRow) {
        alert('Could not find row data');
        return;
    }

    // Check if this is a trio row
    const isTrio = currentRow.Sample.endsWith('-trio');
    const bamPaths = isTrio ?
        (Array.isArray(currentRow.BAM) ? currentRow.BAM : currentRow.BAM.split(',')) :
        [currentRow.BAM];
    const baiPaths = isTrio ?
        (Array.isArray(currentRow.BAI) ? currentRow.BAI : currentRow.BAI.split(',')) :
        [currentRow.BAI];

    if (isLocalFile) {
        // For local files, prompt user to select project directory
        const input = document.createElement('input');
        input.type = 'file';
        input.webkitdirectory = true;
        input.directory = true;
        input.onchange = (e) => {
            const files = Array.from(e.target.files);
            const projectDirectory = files[0].webkitRelativePath.split('/')[0];

            // Find all the files by matching their relative paths
            const bamFiles = bamPaths.map(bamPath =>
                files.find(f => f.webkitRelativePath.replace(`${projectDirectory}/`, '') === bamPath)
            );
            const baiFiles = baiPaths.map(baiPath =>
                files.find(f => f.webkitRelativePath.replace(`${projectDirectory}/`, '') === baiPath)
            );

            if (bamFiles.every(f => f) && baiFiles.every(f => f)) {
                // Read all files using FileReader
                const readers = [...bamFiles, ...baiFiles].map(file => {
                    const reader = new FileReader();
                    return new Promise((resolve) => {
                        reader.onload = () => resolve(reader.result);
                        reader.readAsDataURL(file);
                    });
                });

                Promise.all(readers).then(results => {
                    const bamDataUrls = results.slice(0, bamFiles.length);
                    const baiDataUrls = results.slice(bamFiles.length);
                    initializeIGVViewer(chrom, start, end, bamDataUrls, baiDataUrls, isTrio);
                });
            } else {
                alert(`Could not find the correct BAM/BAI files in the selected directory.\nExpected paths:\n${bamPaths.join('\n')}\n${baiPaths.join('\n')}`);
            }
        };
        input.click();
    } else {
        // For server files, create URLs with the BAM paths
        const baseUrl = window.location.origin + window.location.pathname.replace(/\/[^\/]*$/, '');
        const bamUrls = bamPaths.map(path => baseUrl + '/' + path.replace(/^\//, ''));
        const baiUrls = baiPaths.map(path => baseUrl + '/' + path.replace(/^\//, ''));
        initializeIGVViewer(chrom, start, end, bamUrls, baiUrls, isTrio);
    }
}

// Function to initialize IGV viewer
async function initializeIGVViewer(chrom, start, end, bamUrls, baiUrls, isTrio) {
    // Remove any existing event listeners
    $('#igv-modal').off('shown.bs.modal hidden.bs.modal');
    const container = document.getElementById('igv-viewer');

    // Show the modal
    $('#igv-modal').modal('show');

    // Wait for modal to be fully shown
    $('#igv-modal').on('shown.bs.modal', async function () {
        try {
            // If we don't have a browser instance yet, create one
            if (!window.igvBrowser) {
                const options = {
                    genome: "hg38",
                    locus: `${chrom}:${start}-${end}`,
                    tracks: []
                };
                console.log("Creating IGV browser");
                window.igvBrowser = await igv.createBrowser(container, options);
            } else {
                // If we have an existing browser, just update the locus
                window.igvBrowser.search(`${chrom}:${start}-${end}`);
            }

            // Configure the BAM tracks
            const trackConfigs = bamUrls.map((bamUrl, index) => ({
                name: isTrio ? ["Paternal", "Maternal", "Sample"][index] : "BAM",
                url: bamUrl,
                type: "alignment",
                format: "bam",
                visibilityWindow: 300000000,
                autoHeight: true,
                indexURL: baiUrls[index],
                indexed: true,
                supportsWholeGenome: false,
                groupBy: "tag:HP",
                colorBy: "tag:YC",
                displayMode: "SQUISHED",
                showSoftClips: true,
            }));

            // Load all tracks
            for (const config of trackConfigs) {
                await window.igvBrowser.loadTrack(config);
            }
        } catch (error) {
            console.error("Error with IGV browser:", error);
            alert("Error loading IGV viewer. Please try again.");
        }
    });

    // Clean up when modal is hidden
    $('#igv-modal').on('hidden.bs.modal', function () {
        if (window.igvBrowser) {
            try {
                // Remove all tracks but keep the browser instance
                window.igvBrowser.removeAllTracks();

                // Clean up data URLs if they exist
                if (window.igvDataUrls) {
                    [...window.igvDataUrls.bams, ...window.igvDataUrls.bais].forEach(url => {
                        if (url) URL.revokeObjectURL(url);
                    });
                    window.igvDataUrls = null;
                }
            } catch (error) {
                console.error("Error cleaning up IGV tracks:", error);
            }
        }
    });
}

// Add the viewTrio function
function viewTrio(sample, paternalId, maternalId, event) {
    event.preventDefault();
    event.stopPropagation();

    // Clear any existing filters
    while ($.fn.dataTable.ext.search.length > 0) {
        $.fn.dataTable.ext.search.pop();
    }

    // Add new search function for Sample column (index 5)
    $.fn.dataTable.ext.search.push((settings, data) => {
        const sampleValue = data[5].toLowerCase();  // Sample is column 5
        var sampleName = sample;
        if (sampleName.endsWith("-trio")) {
            sampleName = sampleName.slice(0, -5);
        }
        const searchTerms = [sampleName, sampleName + "-trio", paternalId, maternalId].map(s => s.toLowerCase());
        return searchTerms.some(term => sampleValue === term);
    });

    variant_table.draw();
}

// Add the toggleSpecialInfoFilter function
function toggleSpecialInfoFilter(checkbox) {
    specialInfoFilterActive = checkbox.checked;  // Store the state

    // Find and remove any existing special info filter
    const filterIndex = $.fn.dataTable.ext.search.findIndex(filter =>
        filter.toString().includes('SpecialInfo')
    );
    if (filterIndex !== -1) {
        $.fn.dataTable.ext.search.splice(filterIndex, 1);
    }

    if (specialInfoFilterActive) {
        // Add filter for rows with special info
        $.fn.dataTable.ext.search.push((settings, data, dataIndex) => {
            const rowData = variant_table.row(dataIndex).data();
            return rowData.SpecialInfo !== undefined && rowData.SpecialInfo !== null && rowData.SpecialInfo !== '';
        });
    }

    variant_table.draw();
}

// Add the toggleShowHelp function
function toggleShowHelp(checkbox) {
    showHelpActive = checkbox.checked;
    const helpDiv = document.getElementById('help-div');
    if (helpDiv) {
        helpDiv.style.display = showHelpActive ? 'block' : 'none';
    }
}