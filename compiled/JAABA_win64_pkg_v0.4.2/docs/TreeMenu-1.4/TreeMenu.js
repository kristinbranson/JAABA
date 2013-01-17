/*
 * TreeMenu v1.4
 *
 * Copyright (c) 2006 Mackley F. Pexton.  All rights reserved.
 *
 * This is free software for individual, educational, and non-profit
 * use provided that this copyright notice appears on all copies.
 * Instructions and source code are available at http://www.acmebase.org/tree_menu.
 * Optimized version is for sale at https://order.acmebase.com.
 * Commercial web sites are required to have an inexpensive license.
 * Send correspondence and feedback to: mack_pexton[at]acmebase.org.
 */

/******************************************************************************

TreeMenu v1.4 -- Make menus out of UL/LI tags that open up when clicked.

Tree menus are menus that pop open when a little graphic placed to the left
of the menu item is clicked -- very similar to Microsoft Explorer expanding and
collapsing folders of files.

Usage:

    make_tree_menu(id);
    -or-
    make_tree_menu(id,omit_symbols,no_save_state,singular,no_setup)

where id is the id attribute of the beginning <UL> tag. There are four
optional arguments. The omit_symbols argument, if set to true, omits
adding symbol tags to the beginning of each menu item. The no_save_state
argument prevents the open/close state of the menus from being saved
in a browser cookie. The singular argument, if true, restricts the open
menus to only one per level. The no_setup flag prevents TreeMenu from
scanning the document to setup the menu. This flag assumes that the server
has constructed a complete menu with appropriate symbols, click event handlers,
and class attributes. The optional arguments can be alternately specified
by setting the configuration variables before calling make_tree_menu().

Menus are created out of <UL> tags and their enclosed <LI> tags. The <LI>
tags can contain other <UL> tags. Symbol objects (e.g. <SPAN> tags) are
(optionally) inserted before each <LI> tag. The class name of the symbol
tag is set to one of the TreeMenu class variables: SymbolClassOpen,
SymbolClassClose, or SymbolClassItem depending upon if the <LI> tag has
enclosed <UL> tags and whether the saved open/closed state of the menu's <UL>
tag is 1 or 0. Menu open/closed states are preserved in a cookie with a name
derived from the menu's first <UL> tag's id.

When the user clicks on a symbol tag, the menu is opened or it is closed if it
is already open. The symbol's CSS class is changed to the SymbolClassOpen
or SymbolClassClose appropriately.

To allow more animation, the <LI> tag's class has the class name in the
variable ClassOpen appended to it. This doesn't change the class name,
it adds a new class name. When the menu is closed, the added class name
is removed and the class name in the variable ClassClose added to it.

The type of a symbol's tag inserted before each <LI> tag item is determined
by the class variable SymbolTag.  It can be set to "span" or "div" or it
can be empty. The default value is "span".  A symbol is NOT inserted at the
<LI> tags' beginning if the SymbolTag is set to null or the empty string.
Also, a symbol is NOT inserted at the beginning of the <LI> tags if the 
configuration variable TreeMenu.OmitSymbols is true (default is false).
This allows the server to insert the symbol tags and have this JavaScript
operate the symbols.

You can use regular <A> tags as buttons to open and close menus by defining
an onclick handler like:

  <a href="javascript:;" onclick="TreeMenu.toggle(this)">Submenu</a>

If the button is outside of the menu structure, pass TreeMenu.toggle()
the id of the menu/submenu you want to toggle instead of "this".

  <a href="javascript:;" onclick="TreeMenu.toggle('submenu_id')">Toggle Submenu</a>
  -or-
  <a href="javascript:;" onclick="TreeMenu.show('submenu_id')">Show Submenu</a>
  -or-
  <a href="javascript:;" onclick="TreeMenu.hide('submenu_id')">Hide Submenu</a>

The TreeMenu.show(ul) and TreeMenu.hide(ul) shows and hides one branch
of the tree menu. The TreeMenu.show_all(ul) and TreeMenu.hide_all(ul)
shows and hides every branch of the tree menu. TreeMenu.save_state(ul) saves
the current open/close menu states -- usefull after TreeMenu.show_all()
or TreeMenu.hide_all(). TreeMenu.reset(ul) resets the menu back to the
original default state. It does that by removing the cookie that saves the
menu state. The page needs to be refreshed to actually display the menu in
its original state.

Multiple menus can be made and each menu can have different settings for
the configuration variables.


Changes from TreeMenu v1.3:

- Added use of IMG tags for SymbolTag. This allows the program to change both
the forground images and the background image.  Added TreeMenu.SymbolSrcOpen,
TreeMenu.SymbolSrcClose and TreeMenu.SymbolSrcItem variables.

- Added TreeMenu.Singular flag to restrict open menus to just one per menu
level, and added fourth optional argument to make_tree_menu() to turn on
the flag for an individual menu.

- Added TreeMenu.save_state() utility to save menu open/close states.

- Added TreeMenu.OmitSymbols to refrain from inserting symbols into the list
but yet operate them if found.

- Added TreeMenu.SetupMenu flag to keep TreeMenu from scanning the
document. This allows the server to build the menu instead of TreeMenu
assembling the menu after the document is loaded. This greatly speeds up the
rendering of large menus -- especially with IE6 where accessing the document
objects is so very slow.

******************************************************************************/


/////// Configuration Variables ///////////////////////////

TreeMenu.SymbolTag = 'span';			// symbol inserted at beginning of <LI> tags
//TreeMenu.SymbolTag = '';			// uncomment to disable insertion of symbols
//TreeMenu.SymbolTag = 'img';			// uncomment to use IMG tags for symbols
//TreeMenu.SymbolSrcItem = '';			// url to assign to IMG src attribute of an item
//TreeMenu.SymbolSrcClose = '';			// url to assign to IMG src attribute upon close
//TreeMenu.SymbolSrcOpen = '';			// url to assign to IMG src attribute upon open

TreeMenu.OmitSymbols = false;			// don't insert symbol but do adjust them

TreeMenu.SymbolClassItem  = 'symbol-item';
TreeMenu.SymbolClassClose = 'symbol-close';
TreeMenu.SymbolClassOpen  = 'symbol-open';

TreeMenu.ClassItem  = 'item';			// class name added to <LI> tag's class
TreeMenu.ClassClose = 'close';			// class name added to <LI> tag's class
TreeMenu.ClassOpen  = 'open';			// class name added to <LI> tag's class
TreeMenu.ClassLast  = 'last';			// added to last <LI> and symbol tags' classes

TreeMenu.CookieSaveStates = true;		// flag to use a cookie to save menu state
TreeMenu.CookieExpire = 90;			// days before cookie saving menu states expires

TreeMenu.SetupMenu = true;			// scan document objects to initialize menu

TreeMenu.Singular = false;			// restrict open menus to only one per level

/////// End of Configuration Variables ///////////////////

function make_tree_menu(id,omit_symbols,no_save_state,singular,no_setup) {
	var m = new TreeMenu(id);
	if (omit_symbols) m.OmitSymbols = true;
	if (no_save_state) m.CookieSaveStates = false;
	if (singular) m.Singular = true;
	if (no_setup) m.SetupMenu = false;
	// Setup menus if we are inserting symbols or restoring menu open/close states.
	if (m.SetupMenu) m.setup_symbols();
	return m;
}

/*
 * TreeMenu
 */

function TreeMenu(ul_id) {			// object constructor

	this.top_ul_id = ul_id;
	this.top_ul = document.getElementById(ul_id);

	this.configure();

	// Register menu
	TreeMenu.menus[ul_id] = this;

	return this;
}

/*
 * TreeMenu Class Variables
 */

TreeMenu.menus = [];				// list of defined menus

/*
 * TreeMenu Class Methods
 */

TreeMenu.toggle = function(e) {
	e = TreeMenu.get_ref(e);
	var m = TreeMenu.menus[TreeMenu.get_top_ul(e).id];
	var li = TreeMenu.get_li(e);
	var ul = li.getElementsByTagName("UL")[0];
	if (ul.style.display == "block") {
		m.hide_menu(ul,li,e);
	}
	else {
		if (m.Singular) m.hide_menus_except(li);
		m.show_menu(ul,li,e);
	}

	m.save_menu_states();
}

TreeMenu.show = function(ul) {
	ul = TreeMenu.get_ref(ul);
	var top_ul = TreeMenu.get_top_ul(ul);
	if (! top_ul) return;
	var m = TreeMenu.menus[top_ul.id];
	var li = TreeMenu.get_li(ul);
	m.show_menu(ul,li);
}

TreeMenu.hide = function(ul) {
	ul = TreeMenu.get_ref(ul);
	var top_ul = TreeMenu.get_top_ul(ul);
	if (! top_ul) return;
	var m = TreeMenu.menus[top_ul.id];
	var li = TreeMenu.get_li(ul);
	m.hide_menu(ul,li);
}

TreeMenu.show_all = function(ul) {
	// Show all menus under ul.
	ul = TreeMenu.get_ref(ul);
	var uls = ul.getElementsByTagName("UL");
	for (i = 0; i < uls.length; i++) {
		TreeMenu.show(uls[i]);
	}
}

TreeMenu.hide_all = function(ul) {
	// Hide all menus under ul.
	ul = TreeMenu.get_ref(ul);
	var uls = ul.getElementsByTagName("UL");
	for (i = 0; i < uls.length; i++) {
		TreeMenu.hide(uls[i]);
	}
}

TreeMenu.save_state = function(ul) {
	// Reset menu to original settings.
	ul = TreeMenu.get_ref(ul);
	var m = TreeMenu.menus[TreeMenu.get_top_ul(ul).id];
	m.save_menu_states();
}

TreeMenu.reset = function(ul) {
	// Reset menu to original settings.
	ul = TreeMenu.get_ref(ul);
	var m = TreeMenu.menus[TreeMenu.get_top_ul(ul).id];
	m.reset_menu_states();
}

// Private methods
TreeMenu.get_ref = function(id) {
        if (typeof id == "string") return document.getElementById(id);
	return id;
}

TreeMenu.get_top_ul = function(e) {
	while (e && (e.nodeName != 'UL' || ! e.id || ! TreeMenu.menus[e.id])) e = e.parentNode;
	return e;
}

TreeMenu.get_li = function(e) {
	while (e && e.nodeName != 'LI') e = e.parentNode;
	return e;
}


/*
 * TreeMenu Object Methods
 */

TreeMenu.prototype.configure = function() {

        // Assign global class settings (capitalized variables) to object settings.

        var v,c;
        for (v in TreeMenu) {
                c = v.substr(0,1);
                if (c == c.toUpperCase()) {
                        this[v] = TreeMenu[v];
                }
        }
}

TreeMenu.prototype.setup_symbols = function() {

	// Insert open/close symbols at the beginning of the menu items
	// and open or close menus like they were previously.

	var states = this.get_menu_states();

	var index = 0;
	var ul, li, symbol, islast = false;
	var ul_elements, li_elements = this.top_ul.getElementsByTagName("LI");
	for(var i=0; i < li_elements.length; i++) {
		li = li_elements[i];

		if (this.ClassLast) islast = this.is_last_item(li);

		ul_elements = li.getElementsByTagName("UL");
		if(ul_elements.length > 0) {
			// Submenus
			if (this.SymbolTag && ! this.OmitSymbols) {
				symbol = document.createElement(this.SymbolTag);
				if (this.ClassLast && islast) symbol.className = this.ClassLast;
				symbol.onclick = function() { TreeMenu.toggle(this); };
				li.insertBefore(symbol, li.firstChild);
			}

			ul = ul_elements[0];
			if (states[index] == '1') this.show_menu(ul,li);
			else                      this.hide_menu(ul,li);
			index++;
		}
		else {
			// Menu item
			if (this.SymbolTag && ! this.OmitSymbols) {
				symbol = document.createElement(this.SymbolTag);
				if (this.SymbolClassItem)
					symbol.className = this.SymbolClassItem;
				if (this.SymbolSrcItem)
					symbol.src = this.SymbolSrcItem;
				if (this.ClassLast && islast)
					symbol.className += ' ' + this.ClassLast;
				li.insertBefore(symbol, li.firstChild);
			}

			if (this.ClassItem) li.className += ' ' + this.ClassItem;
		}

		if (islast) li.className += ' ' + this.ClassLast;
	}
}

TreeMenu.prototype.is_last_item = function(e) {
	// Check if element is the last LI element in the list.
	e = e.nextSibling;
	// Get next element (Mozilla puts text nodes at same level here).
	while (e && e.nodeType != 1) e = e.nextSibling;
	return e ? false : true;
}

TreeMenu.prototype.get_menu_states = function() {
	var cookie = getCookie("tm_" + this.top_ul_id);
	if (cookie) return cookie.split('x');
	return [];
}

TreeMenu.prototype.save_menu_states = function() {

	// Save all menu and submenu open/close states in a cookie

	if (! this.CookieSaveStates) return;

	var states = [];
	var ul_elements, li_elements = this.top_ul.getElementsByTagName("LI");
	for(var i=0; i < li_elements.length; i++) {
		ul_elements = li_elements[i].getElementsByTagName("UL");
		if (ul_elements.length > 0) {
			states[states.length] = ul_elements[0].style.display == "block" ? 1 : 0;
		}
	}

	var expire_date = new Date((new Date().getTime()) + this.CookieExpire*24*60*60*1000);
	setCookie("tm_" + this.top_ul_id, states.join('x'), expire_date, '/');
}

TreeMenu.prototype.reset_menu_states = function() {

	// Reset all menu and submenu open/close states  (delete cookie)

	var expire_date = new Date((new Date().getTime()) - 1000);		// set to past time
	setCookie("tm_" + this.top_ul_id, '', expire_date, '/');
}

TreeMenu.prototype.add_remove_class = function(e,add_class,remove_class) {
	if (e) {
		if (remove_class)
			e.className = e.className.replace(remove_class,'');
		if (add_class && ! e.className.match( (new RegExp("\\b"+add_class+"(\\s.*)?")) ) ) {
			e.className += ' ' + add_class;
		}
	}
}

TreeMenu.prototype.show_menu = function(ul,li,e) {
	ul.style.display = 'block';

	this.add_remove_class(li,this.ClassOpen,this.ClassClose);

	if (this.SymbolTag) {
		var symbol = li.getElementsByTagName(this.SymbolTag)[0];
		this.add_remove_class(symbol,this.SymbolClassOpen,this.SymbolClassClose);
		if (this.SymbolSrcOpen) symbol.src = this.SymbolSrcOpen;
	}

	// Following case is for toggle buttons disassociated with menu structure.
	this.add_remove_class(e,this.SymbolClassOpen,this.SymbolClassClose);
}

TreeMenu.prototype.hide_menu = function(ul,li,e) {
	ul.style.display = 'none';

	this.add_remove_class(li,this.ClassClose,this.ClassOpen);

	if (this.SymbolTag) {
		var symbol = li.getElementsByTagName(this.SymbolTag)[0];
		this.add_remove_class(symbol,this.SymbolClassClose,this.SymbolClassOpen);
		if (this.SymbolSrcClose) symbol.src = this.SymbolSrcClose;
	}

	// Following case is for toggle buttons disassociated with menu structure.
	this.add_remove_class(e,this.SymbolClassClose,this.SymbolClassOpen);
}

TreeMenu.prototype.hide_menus_except = function(li) {
	// Hide other menus at same level as li.
	var n;
	var re = new RegExp('\\b' + this.ClassOpen + '\\b');
	for (var i = 0; i < li.parentNode.childNodes.length; i++) {
		n = li.parentNode.childNodes[i];
		if (n == li || n.nodeType != 1) continue;
		if (n.className.match(re)) this.hide_menu(n.getElementsByTagName("UL")[0],n);
	}
}

/*
 * Classic Cookie functions
 */

function setCookie(name, value, expires, path, domain, secure) {
	document.cookie	= name + "=" + escape(value) +
	  (expires	? "; expires=" + expires.toGMTString()	: "") +
	  (path		? "; path=" + path			: "") +
	  (domain	? "; domain=" + domain			: "") +
	  (secure	? "; secure"				: "");
}

function getCookie(name) {
	var dc = document.cookie;
	var prefix = name + "=";
	var begin = dc.indexOf("; " + prefix);
	if (begin == -1) {
		begin = dc.indexOf(prefix);
		if (begin != 0) return null;
	}
	else {
		begin += 2;
	}
	var end = document.cookie.indexOf(";", begin);
	if (end == -1) end = dc.length;
	return unescape(dc.substring(begin + prefix.length, end));
}
