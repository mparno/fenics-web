<ul class="simple">
  <li>
    <a onclick="javascript:ShowHide('HiddenDivMain')" href="javascript:;">Show/Hide</a>
  </li>
  <div class="mid" id="HiddenDivMain" style="DISPLAY: none" >
    <iframe src="http://www.fenicsproject.org:8080"
	    frameBorder="0"
	    height="300"
	    width="100%"
	    align="center"
	    name="dolfin-main">
      <a href="http://www.fenicsproject.org:8080">http://www.fenicsproject.org:8080</a>
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
