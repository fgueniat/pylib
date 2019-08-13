from bs4 import BeautifulSoup
import urllib2
from selenium import webdriver
import time
import div_tools as divt
import os

#not invisible
#driver = webdriver.Firefox()
#invisible

class wishlist():
    def __init__(self,path=None,name = None,loginpath=None,login=None,password=None):
        

        if name is None:self.name = os.path.expanduser('~')+'/temp/wishlist.data'
        else: self.name = name
        self.gixen_ids = []
        self.results = {}
        
        if loginpath is not None:
            try:
                r=divt.load(loginpath)
            except:
                'pass incorrect'
        else:
            loginpath = os.path.expanduser('~')+'/temp/login.data'
            try:
                r=divt.load(loginpath)
            except:pass

        if login is not None: 
            self.login = login
        else:
            try:self.login = r['login']
            except: self.login=raw_input('login needed: ')

        if password is not None: 
            self.password = password
        else:
            try:self.password = r['password']
            except: self.password=raw_input('password needed: ')

        object_list = {
                ## Black
                    'liliana':
                    {   'active':True,
                            'url':'http://www.ebay.com/sch/i.html?_from=R40&_sacat=0&_nkw=liliana+of+the+veil&_sop=15',
                            'banned_words':['altered','bag','playmat'],'forced_words':None,
                            'max':50},            
        #
                    'obliterator':
                    {   'active':True,
                            'url':'http://www.ebay.com/sch/i.html?_odkw=sheoldred+whispering+one&_sop=15&_osacat=0&_from=R40&_trksid=p2045573.m570.l1311.R1.TR11.TRC1.A0.H0.Xphyrexi.TRS0&_nkw=phyrexian+obliterator&_sacat=0','banned_words':None,'forced_words':None,
                            'max':15},
                ## Colorless
                    'forcefield':
                        {
                            'active':True,
                            'url':'http://www.ebay.com/sch/i.html?_odkw=forcefield+mtg+-collectors+-international+-+Collector&_sop=15&_osacat=0&_from=R40&_trksid=p2045573.m570.l1313.TR0.TRC0.H0.Xforcefield+unlimitedmtg+.TRS2&_nkw=forcefield+unlimited+mtg+&_sacat=0',
                            'max':120,
                            'banned_words':['international','collector','ce'],
                            'forced_words':['forcefield'],},
        #
                    'lion eye diamond':
                        {
                            'active':True,
                            'url':'http://www.ebay.com/sch/i.html?_from=R40&_trksid=p2380057.m570.l1311.R2.TR4.TRC0.A0.H1.Xlion+eye.TRS0&_nkw=lion+eye+diamond+mtg&_sacat=0','banned_words':None,'forced_words':None,
                            'max':70},
                ## Green
                    'life from the loam':
                        {
                            'active':False,
                            'url':'http://www.ebay.com/sch/i.html?_from=R40&_trksid=p2380057.m570.l1313.TR0.TRC0.H0.TRS0&_nkw=Life+from+the+Loam+&_sacat=0','banned_words':None,'forced_words':None,
                            'max':10},            
                ## Blue
                        'jace':
                         {
                             'active':True,
                             'url':'http://www.ebay.com/sch/i.html?_from=R40&_trksid=p2380057.m570.l1311.R1.TR11.TRC1.A0.H0.Xjace+the+m.TRS0&_nkw=jace+the+mind+sculptor&_sacat=0','banned_words':None,'forced_words':None,
                            'max':45},

                    }

        loadsuccess = False
        if path is not None:
            try:
                self.load(path)
                loadsuccess = True
            except:
                print 'path is incorrect. Regular path is tried'

        if loadsuccess is False:
            try: 
                self.load()
                loadsuccess = True
            except: 
                print 'regular path failed, hard coded list used'
        if loadsuccess is False:
            self.list   = {}
            for key in object_list.keys():
                self.add(name = key,price = object_list[key]['max'],url =object_list[key]['url'],forced = object_list[key]['forced_words'],active = object_list[key]['active'],banned = object_list[key]['banned_words'])


    def load(self,path=None):
        if path is None: path = self.name
        object_list = divt.load(path)
        self.list = {}
        for key in object_list.keys():
            self.add(name = key,price = object_list[key]['max'],url =object_list[key]['url'],forced = object_list[key]['forced_words'],active = object_list[key]['active'],banned = object_list[key]['banned_words'])
 
    def delete(self,card):
        try: 
            del self.list[card]
            self.save()
        except: print 'item not present'
            
    def save(self,name=None):
        if name is None:name=self.name
        divt.save(obj=self.list,filename = name,verbose=False)

    def set(self,url=None,name=None,price=None,banned=None,forced=None,active=None):
        if name not in self.cards:
            print 'card not present'
            print 'cards are: '
            print self.cards
            return -1

        if url is not None: self.list[name]['url'] = url
        if price is not None: self.list[name]['price'] = price
        if banned is not None: self.list[name]['banned_words'] = banned
        if forced is not None: self.list[name]['forced_words'] = forced
        if active is not None: self.list[name]['active'] = active

        self.save()

       
    def add(self,url=None,name=None,price=None,banned=None,forced=None,active=True):
        if url is None: url=raw_input('url of the search: ')
        if price is None: price=raw_input('max price of the search: ')
        if name is None:name=raw_input('name of the search: ')
        if name in self.list.keys():
            print name + 'is updated'
        if name is not None:
            self.list[name] = {
                'active':active,
                'banned_words':banned,
                'forced_words':forced,
                'max':price,
                'url':url,
                }
            self.save()
        else: print 'no name - no save'

    def isactive(self,item=None):
        if item is None: return False
        try:
            return self.list[item]['active']
        except:
            print 'error in isactive (with key?)'
            return False
    def price(self,item=None):
        if item is None: return 10000000000
        try:
            return float(self.list[item]['max'])
        except:
            print 'error in price (with key?)'
            return 1000000
    def url(self,item=None):
        if item is None: return ''
        try:
            return self.list[item]['url']
        except:
            print 'error in url (with key?)'
            return '' 
    def forced_words(self,item=None):
        if item is None: return None
        try:
            return self.list[item]['forced_words']
        except:
            print 'error in forced_words (with key?)'
            return None
    
    def banned_words(self,item=None):
        if item is None: return None
        try:
            return self.list[item]['banned_words']
        except:
            print 'error in banned_words (with key?)'
            return None

    def update_gixen(self):
        #define the driver: phantom for invisible, if note Firefox.
        driver = webdriver.PhantomJS()
        #open url
        driver.get("https://www.gixen.com/index.php")
        #fill forms
        element = driver.find_element_by_name("username")
        element.send_keys(self.login)
        element = driver.find_element_by_name("password")
        element.send_keys(self.password)
        element = driver.find_element_by_name('Submit')
        element.click()
        time.sleep(5)
        #read the source for beautifulsoup
        html = driver.page_source
        driver.close()
        #Hence the soup
        self.soup_gixen = BeautifulSoup(html,'lxml')
        #found the items
        soup_ids = self.soup_gixen.find_all('input',type='hidden',id='edititemid')

        self.gixen_ids = []
        #add the prices
        for id_ in soup_ids:
            value = id_['value']
            try:
                value = int(value)
                self.gixen_ids.append(value)
            except:
                'error in add gix value'
        if len(self.gixen_ids)<1:print 'either no items or login failed'
        print "gixen updated"

    def update(self): 
        self.update_gixen()
        self.scrap_ebay()
    
    def scrap_(self,name=None):self.scrap_ebay(name=name)

    def scrap_ebay(self,name=None):
        
        if name is not None:
            if name not in self.list.keys():
                print 'add card to the list first'
                return -1
            else: list_of_cards = [name]
        else: list_of_cards = self.cards
           
        self.results = {}
        for card in list_of_cards:
            if self.isactive(card) is True:
                url2load = urllib2.urlopen(self.url(card))
                read_html = url2load.read()
                url2load.close()

                soup = BeautifulSoup(read_html,'lxml')
                # check prices
                indices = []
                prices_temp = []
                soup_prices = soup.find_all('li',class_='lvprice prc')
                for ind_,price_raw in enumerate(soup_prices):
                    #clean the prices
                    try:
                        price = str(price_raw.get_text())
                        price = price.replace('\n','')
                        price = price.replace('$','')
                        price = price.replace(',','')
                        price = float(price)
                    except:
                        print 'error'
                        print price
                        price = 10000000.
                        
                    if price < self.price(card):
                        indices.append(ind_)
                        prices_temp.append(price)

                # check data
                ebay_ids = []
                urls = []
                titles = []
                prices = []
                soup_ids = soup.find_all('div',class_ = 'lvpic pic img left')
                soup_urls = soup.find_all('a',class_ = 'vip')

                try:
                    banned_words = ['repack', 'proxy', 'pack','custom']
                    forced_words = [card]
                    if self.banned_words(card) is not None: banned_words = banned_words + self.banned_words(card) 
                    if self.forced_words(card) is not None: forced_words = forced_words + self.forced_words(card)
                except:
                    print 'error:'
                    forced_word = []
                    print card
                    print [card]


                for i,ind_ in enumerate(indices):
                    try:
                        # parse data
                        title = soup_urls[ind_]['title'].replace('Click this link to access ','')
                        url_ = soup_urls[ind_]['href']
                        id_ = int(soup_ids[ind_]['iid'])
                        
                        isvalid = True
                        #check titles
                        for word in forced_words:
                            if word.lower() not in title.lower():
                                isvalid = False
                                #print title
                        for word in banned_words:
                            if word.lower() in title.lower():
                                isvalid = False
                                #print title
                        ############# TO DO : check url ?
                        if id_ in self.gixen_ids:
                                isvalid = False
                        if isvalid is True: #test is passed
                            #title
                            titles.append(title)
                            #url
                            url_ = soup_urls[ind_]['href']
                            urls.append(url_)
                            #ID
                            id_ = int(soup_ids[ind_]['iid'])
                            ebay_ids.append(id_)
                            prices.append(prices_temp[i])
                    except:
                        error = 'an error happens. Probably len(prices)>len(data)'
                result = {'urls':urls,'ids':ebay_ids,'prices':prices}
                self.results[card] = result
            self.show(card)
        
    def show(self,name=None):
        if name is not None:
            if name not in self.list.keys():
                print 'add card to the list first'
                return -1
            else: list_of_cards = [name]
        else: list_of_cards = self.cards

        for card in list_of_cards:
            if self.isactive(card) is True:
                print '#########################################'
                print card
                for i in xrange(len(self.results[card]['prices'])):
                    print card
                    print '## prix: ' + str(self.results[card]['prices'][i]) + ' et url:'
                    print self.results[card]['urls'][i]
                    print '## id: ' + str(self.results[card]['ids'][i])

    @property
    def cards(self):
		return self.list.keys()

if __name__ == '__main__':
    print "This part is executed only when %s is executed rather than imported" % __file__
    liste = wishlist()
    liste.update()

