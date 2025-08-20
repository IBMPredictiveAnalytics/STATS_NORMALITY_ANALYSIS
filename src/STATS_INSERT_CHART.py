# anonymous extension command to insert a list of charts into the Viewer

# history
# 13-Nov-2024  Original vesion
# 09-Dec-2024  update to handle multiple charts on one call via a file of names
# 27-apr-2025  minor edits

# Author: Jon K. Peck

# No help file is provided for this extension

# USAGE:
# STATS INSERT CHART = filespec POSITION=where to insert HEADER="header text" OUTLINELABEL="text for outline"
# HIDELOG = YES or NO.
 
import spss, SpssClient
from extension import Template, Syntax, processcmd

def doinsertcharts(chartlist=None, header=None, outlinelabel=None, labelparm=None, hidelog=False):
    """Insert one or more charts into the Viewer
    
    chartlist is the filespec of the file holding list of chart to insert
    header is the outline title for the charts
    outline label is the label for the item
    labelparm is a list of appends to the outline label or None
    hidelog specifies that the most recent log block in the Viewer should be closed
    """
    # debugging
            # makes debug apply only to the current thread
    #try:
        #import wingdbstub
        #import threading
        #wingdbstub.Ensure()
        ####wingdbstub.debugger.SetDebugThreads({threading.get_ident(): 1})
    #except:
        #pass

    if chartlist is not None:
        try:
            with open(chartlist) as f:
                charts = f.readlines()
        except:
            print("CHARTLIST file not found")
            raise
        
        charts = [item[:-1] for item in charts]  # eliminate newlines
    else:
        charts = None

    SpssClient.StartClient()
    if chartlist is not None and (header is None or outlinelabel is None):
        raise ValueError("Missing keyword value")
    doc = SpssClient.GetDesignatedOutputDoc()
    itemlist = doc.GetOutputItems()
    # Get the root header item
    root = itemlist.GetItemAt(0).GetSpecificType()
    theHeader = doc.CreateHeaderItem(header)
    root.InsertChildItem(theHeader, root.GetChildCount())
    headerItem = root.GetChildItem(root.GetChildCount()-1)
    if not headerItem.GetType() ==  SpssClient.OutputItemType.HEAD:
        root.RemoveChildItem(root.GetChildCount()-1)
        headerItem =  root.GetChildItem(root.GetChildCount()-1)
    headerItem = root.GetChildItem(root.GetChildCount()-1).GetSpecificType()

    for position, chart in enumerate(charts):
        # Create a new chart item
        if labelparm:
            lbl = outlinelabel + str(labelparm[position])
        else:
            lbl = outlinelabel
        outitem = doc.CreateImageChartItem(chart,f"{lbl}")
        # Append the new item to the header item
        headerItem.InsertChildItem(outitem, position)
        
    if hidelog:
        hidethelog(itemlist)
    SpssClient.StopClient()
    
def hidethelog(itemlist):
    """Hide the most recent log block in the Viewer (if any)

    itemlist is the list of items currently in the Viewer"""
    
    itemkt = itemlist.Size()
    for i in range(itemkt-1, 0, -1):
                    item = itemlist.GetItemAt(i)
                    if item.GetType() == SpssClient.OutputItemType.LOG:
                        item.SetVisible(False)
                        break
    
def Run(args):
    """Execute the STATS INSERT CHART command"""
    
    args = args[list(args.keys())[0]]
    oobj = Syntax([
        Template("CHARTLIST", subc="",  ktype="str", var="chartlist", islist=False),
        Template("LABELPARM", subc="",  ktype="int", var="labelparm", islist=True),
        Template("HEADER", subc="", ktype="literal", var="header", islist=False), 
        Template("OUTLINELABEL", subc="", ktype="literal", var="outlinelabel", islist=False),
        Template("HIDELOG", subc="", ktype="bool", var="hidelog", islist=False)
    ])
        
    #enable localization
    global _
    try:
        _("---")
    except:
        def _(msg):
            return msg

    # A HELP subcommand overrides all else
    if "HELP" in args:
        #print helptext
        helper()
    else:
        processcmd(oobj, args, doinsertcharts)

def helper():
    """open html help in default browser window
    
    The location is computed from the current module name"""
    
    import webbrowser, os.path
    
    path = os.path.splitext(__file__)[0]
    helpspec = "file://" + path + os.path.sep + \
         "markdown.html"
    
    # webbrowser.open seems not to work well
    browser = webbrowser.get()
    if not browser.open_new(helpspec):
        print(("Help file not found:" + helpspec))
try:    #override
    from extension import helper
except:
    pass        
