<html>
  <head>
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load('current', {'packages':['sankey']});
      google.charts.setOnLoadCallback(drawChart);

      function drawChart() {
        var data = new google.visualization.DataTable();
        data.addColumn('string', 'from');
        data.addColumn('string', 'to');
        data.addColumn('number', 'cells');
        data.addRows([
          {{ matrix }}
        ]);

        // Sets chart options.
        var options = {
          width: 1200,
          hieght: 1200,
          sankey: {
            node: {
                    colors: [{{ colors }}]
                },
            link: {
                    colorMode: 'gradient'
                }
            }
        };

        // Instantiates and draws our chart, passing in some options.
        var chart = new google.visualization.Sankey(document.getElementById('sankey_basic'));
        chart.draw(data, options);
      }
    </script>
  </head>
  <body>
    <div id="sankey_basic" style="width: 6000px; height: 3000px;"></div>
  </body>
</html>