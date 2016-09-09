
node["data"]["id"]  # event ID
node["data"]["date"] # the date of the event
node["data"]["title"] # the title of the event (letters in lower case)
node["data"]["catg"] # a list of categories 
node["data"]["concept"] # a list of concepts
node["data"]["summary"] # the summary of the event
node["data"]["location"] # the location this event occurs
node["data"]["country"] # the country this event occurs
node["data"]["score"] # the social score (hotness) of the event

self.G.graph["conceptDict"] # dictionary of concepts. Key: concept; Value: a list of tuples (event ID, date)
self.G.graph["catgDict"] # dictionary of categories. Key category; Value: a list of categories.  
# EP_prediction
