window.onload = function() {
    $.ajax({
	url: 'https://www.googleapis.com/plus/v1/people/105550716956576029273/activities/public?key=INSERT_KEY_HERE',
	success: function(news) {
	    $('#news').append('<p>' + news.items[0].object.content + '</p>');
	    $('#news').append('<small>' + news.items[0].actor.displayName + '&mdash;' + news.items[0].published + '</small>');
	},
	dataType: 'json'
    });
}
