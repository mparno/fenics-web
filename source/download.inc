<script type="text/javascript">
<!--
BrowserDetect.init();

document.write('<p>You appear to be using\
<li> browser: ' + BrowserDetect.browser + '</li>\
<li> version: ' + BrowserDetect.version + '</li>\
<li> os:      ' + BrowserDetect.OS + '</li>\
<li> distro:  ' + BrowserDetect.Distro + '</li>\
</ul>\
</p>');

document.write('<p>Perhaps I could interest you in the <a href="foo-' + BrowserDetect.Distro + '.tar.gz">download page for ' + BrowserDetect.OS + ' ' + BrowserDetect.Distro + '</a>.</p>');
distro = BrowserDetect.Distro;
if (distro == "Hardy"    ||
    distro == "Intrepid" ||
    distro == "Jaunty"   ||
    distro == "Karmic"   ||
    distro == "Lucid"    ||
    distro == "Maverick" ||
    distro == "Natty")
{
document.write('<p>Hooray! You\'re using Ubuntu ' + distro + '. You are eligible for <a href="foo-' + distro + '-packages.deb">the packages</a>.</p>');
}
-->
</script>
