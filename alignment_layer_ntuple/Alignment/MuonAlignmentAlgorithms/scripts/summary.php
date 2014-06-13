<html>
<body>
<?php

$dirname = ".";
$filename = "summary.php";
$indexname = "_index.php";

echo "<html>\n";
$html = "";


$objects = new RecursiveIteratorIterator(new RecursiveDirectoryIterator($dirname), RecursiveIteratorIterator::SELF_FIRST);
$plotType = str_replace(".png","",$_GET["plot"]);
$chambers = array();
foreach($objects as $file) {
    if($file->isDir() ) {
        $depth = substr_count($file, '/');
        if($depth == 4) {
            $chambers[] = $file; // weird syntax to APPEND
            $curimg = $file . "/" . $plotType . ".png";
            $html .= "\t<a href='$curimg'><img src='$curimg' height=45% /></a>\n";
        }
    }
}

echo "<table width='100%' border=0 cellspacing=0 cellpadding=5>\n"; // begin main table
//echo "<table border=1 style='display: inline-block;'>\n";
//echo "<table border=1 width: 50%; float: center;'>\n";
echo "<tr>\n";
echo "<td width='50%'>\n";
echo "<table width='100%' border=0 cellspacing=0 cellpadding=1>\n"; // begin 1st table
echo "<tr><td><b>Chambers</b> (takes you to chamber's index page)</td></tr>\n";
foreach($chambers as $chamber) {
    echo "<tr>\n";
    echo "<td>";
    echo "<a href='$chamber/$indexname'>$chamber</a>";
    echo "</td>\n";
    echo "</tr>\n";
}
echo "</table>\n"; // end 1st table
echo "</td>\n";

echo "<td width='50%'>\n";
echo "<table width='100%' border=0 cellspacing=0 cellpadding=1>\n"; // begin 2nd table
echo "<tr><td><b>Plot types</b> (displays 1 plot type for all chambers here)</td></tr>\n";
$plots = scandir($chambers[0]);
foreach($plots as $plot) {
    if(  strpos($plot,".png") &&
        !strpos($plot,"LAY2") &&
        !strpos($plot,"LAY3") &&
        !strpos($plot,"LAY4") &&
        !strpos($plot,"LAY5") &&
        !strpos($plot,"LAY6")  )
    {
        echo "<tr>\n";
        echo "<td>";
        echo "\t<a href='$filename?plot=$plot'>$plot</a><br>\n";
        echo "</td>\n";
        echo "</tr>\n";
    }

}
echo "</table>\n"; // end 2nd table
echo "</td>\n";
echo "</tr>\n";


echo "</table>\n"; // end main table

echo $html;
echo "</html>";

?>
</body>
</html>
