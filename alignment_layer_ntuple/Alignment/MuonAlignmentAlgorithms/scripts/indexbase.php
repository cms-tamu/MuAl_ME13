<html>
<body>
<?php

$dirname = ".";

$files = scandir($dirname);

echo "<center><br>\n";
foreach($files as $text) {
    if(strpos($text, '.txt')) {
        echo str_replace(".txt","",$text) . ": " . file_get_contents($text) . "<br>\n";
    }
}

echo "\n</p>\n";

foreach($files as $curimg){
    if (strpos($curimg,'.png') || strpos($curimg,'.jpg')) {
        echo "\t<a href='$curimg'><img src='$curimg' height=45% /></a>\n";
    }
}
echo "</center>\n";
?>
</body>
</html>
