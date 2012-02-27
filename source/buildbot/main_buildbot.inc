<ul class="simple">
  <li>
    <a onclick="javascript:ShowHide('HiddenDivTrunk')" href="javascript:;">Show/Hide</a>
    (<a href="https://code.launchpad.net/~dolfin-core/dolfin/trunk">lp:dolfin/trunk</a>)
  </li>
  <div class="mid" id="HiddenDivTrunk" style="DISPLAY: none" >
    <iframe src="http://fenicsproject.org:8090"
	    frameBorder="0"
	    height="400"
	    width="100%"
	    align="center"
	    name="dolfin-trunk">
      <a href="http://fenicsproject.org:8090">http://fenicsproject.org:8090</a>
    </iframe> 
  </div>
</ul>

<ul class="simple">
  <li>
    <a onclick="javascript:ShowHide('HiddenDiv1.0.x')" href="javascript:;">Show/Hide</a>
    (<a href="https://code.launchpad.net/~dolfin-core/dolfin/1.0.x">lp:dolfin/1.0.x</a>)
  </li>
  <div class="mid" id="HiddenDiv1.0.x" style="DISPLAY: none" >
    <iframe src="http://fenicsproject.org:8080"
	    frameBorder="0"
	    height="400"
	    width="100%"
	    align="center"
	    name="dolfin-1.0.x">
      <a href="http://fenicsproject.org:8080">http://fenicsproject.org:8080</a>
    </iframe> 
  </div>
</ul>

<script language="JavaScript">
function ShowHide(divId)
{
    if(document.getElementById(divId).style.display == 'none')
    {
        document.getElementById(divId).style.display='block';
    }
    else
    {
        document.getElementById(divId).style.display = 'none';
    }
}
</script>
