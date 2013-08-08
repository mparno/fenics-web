window.onload = function() {
    $.ajax({
	url: 'https://www.googleapis.com/plus/v1/people/105550716956576029273/activities/public?key=AIzaSyDmcx1Nx8NOiFL_sSwfgq6ovjfJRl6Kf8g',
	success: function(news) {
            for (var i = 0; i < 7; i++)
            {
	        $('#news').append('<p>');
		$('#news').append('<a href="' + news.items[i].url + '"><img class="avatar" src="' + news.items[i].actor.image.url + '"></a>');
	        $('#news').append(news.items[i].title +
                                  ' ' +
                                  '<a href="' + news.items[i].url + '"><br>Read more &raquo;</a>');
	        $('#news').append('<small><br>' +
                                  'Posted by <a href="' + news.items[i].actor.url + '">' + news.items[i].actor.displayName + '</a>' +
                                  '&mdash;' + news.items[i].published +
                                  '</small>');
	        $('#news').append('</p>');
            }
	},
	dataType: 'json'
    });
}
