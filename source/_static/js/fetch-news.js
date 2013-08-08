window.onload = function() {
    $.ajax({
	url: 'https://www.googleapis.com/plus/v1/people/105550716956576029273/activities/public?key=AIzaSyDmcx1Nx8NOiFL_sSwfgq6ovjfJRl6Kf8g',
	success: function(news) {
            for (var i = 0; i < 7; i++)
            {
	        $('#news').append('<p>');
		$('#news').append('<a href="' + news.items[i].url + '"><img class="avatar" src="' + news.items[i].actor.image.url + '"></a>');
	        $('#news').append(news.items[i].title);
		var pubdate = new Date(news.items[i].published);
	        $('#news').append('<small><br>' +
                                  '<a href="' + news.items[i].url + '">Read more</a> &mdash; Posted by <a href="' + news.items[i].actor.url + '">' + news.items[i].actor.displayName + '</a>' +
                                  ' on ' + pubdate.getFullYear() + '-' + pubdate.getMonth() + '-' + pubdate.getDate() +
                                  '</small>');
	        $('#news').append('</p>');
            }
	},
	dataType: 'json'
    });
}
