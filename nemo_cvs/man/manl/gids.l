.TH GIDS 1NEMO "22 July 1994"
.SH NAME
gids \- Groningen Image Display Server
.SH SYNOPSIS
\fBgids\fP [-option value]...
.SH DESCRIPTION
\fBgids\fP is an X-windows display server program  
that can display and manipulate images, in particular store an unlimited
amount of images from a cube and display them in \fImovie\fP and 
\fIslice\fP mode. (see also \fIds9(lNEMO)\fP for a simple standalone
version now in more common use within NEMO)
.PP
Clients talk to this server using the \fIgdi\fP library. The server can
also send information back to the client, and it is up to the client to
interpret this information. The server also understand the PGPLOT
plotting routines.
.SH OPTIONS
Command line options are optional and of the form \fI-option value\fP.
They are:
.TP 20
\fB-display\fP \fIserver\fP
Display server to run on. Usualy X-windows syntax. By default, the
\fBDISPLAY\fP environment variable is used.
.TP
\fB-memory\fP \fIsize\fP
Size, in kilo-bytes, to allocate for image frames. These are used to
store images in movie and slice mode. The maximum number of images is
displayed in the \fBREC=\fP (third) window on the GIDS server. By default,
space for 3 512*512 images is allocated (\fIsize=768\fP).
.TP
\fB-file\fP \fIsetup_file\fP
Setup file, by default \fB~/.default_display\fP, which stores information
about the server.
.TP
\fB-keyboard\fP \fIkeyboard\fP
Also controlled by the environment variable \fBALT_KEYBOARD\fP, 
if present.
.SH RESOURCES
\fBgids\fP allows the user to specify some defaults in the \fB.Xdefaults\fP
or \fB.Xresources\fP file. The parameters which can be set are the
following:
.nf
.ta +0.5i +2i
gids.calibration_font:  	font for calibration window.
gids.help_file:            	path to help file.
gids.help_font:            	font for help window.
gids.info_font:            	font for info windows.
gids.menu_font:            	font for menu windows.
gids.display_height:       	default height of display window.
gids.display_width:        	default width of display window.
gids.start_lut:            	color lookup table at program start
                    		 (use name of LUT as from color menu).
gids.scaling:              	color transfer funtion as program
                		 start (use name of scaling funtion as from scale menu).
gids.inverse_scaling:      	True for negative slope of color transfer function.
gids.autoraise:            	Raise window when new image or plot is arriving [False].

.SH MENUS
Here is a detailed listing of the online-available help menu. It is activated
in context-sensitive mode by pressing the middle mouse button on any menu
option.
.nf
(0)           HELP

              GIDS is the Groningen Image Display Server. It will display
              two dimensional images which have to be loaded into GIDS by a
              special program (i.e. the GIPSY program VIEW).

              MOUSE BUTTONS:
              Here we assume that button1 is the left-most mouse button,
              button2 the middle mouse button and button3 the right-most
              mouse button.

              To remove the HELP window, press HELP again with button1 or
              button3 or press a mouse button in the help window.

              DISPLAY WINDOW:
              In the display window the images are shown. At startup the
              display window is the right-most window with a size of 512 x 512.
              It is possible to resize this window. The minimum size is
              128 x 128, the maximum size depends on the size of the screen.
              Each resizing of the display window will remove all previous
              recordings.

              WEDGE WINDOW:
              To the left of the display window the wedge window is placed.
              This window displays the current color lookup table.

              CALIBRATION WINDOW:
              To the left of the wedge window the calibration window is placed.
              Here the data values as a function of current color are
              displayed graphically.

              INFO WINDOWS:
              GIDS has 7 info windows which are placed in the upper-left
              window of GIDS (white on black). From top to bottom, the first
              two display info related to the image currently in the display
              window. This information depends on what the program which loaded
              the image into GIDS found necessary to show the user.
              The third window displays information about the recordings
              (REC= MR, NR, CR), where MR is the maximum number of recordings
              the server can handle, NR the number of images recorded and
              CR the current record number in the display window. Record 0 is
              special since here always the last image loaded into GIDS is
              stored. It is only possible to (de)zoom this image. Pressing
              the ETCETERA menu and then the RECORDING menu will show a
              list of recording options.
              The fourth window displays information about the graphics planes
              (MASK=  M ON=  N), where M is the current mask of graphics planes
              on which the graphics operation is valid and N the mask for the
              planes which are displayed. GIDS has 4 graphics planes, which are
              numbered 1, 2, 4 and 8. If more than one graphics plane is wanted,
              simply add the numbers. Pressing the ETCETERA menu and then the
              REGION menu will give a list of graphics plane operations.
              The fifth window displays the data value at the position of the
              cursor in the display window (VAL=   N).
              The sixth window will display the grids at the position of the
              cursor in the display window (POS=   X,   Y).
              The seventh window is an input window which displays the text
              entered via the keyboard. Some menu actions need the text in
              this window.

              MENU WINDOWS:
              Below the info windows menu windows are placed.
              Moving the pointer to a menu window will 'highlight' it, clicking
              mouse button1 or button3 will activate it. Clicking button2
              will give help on the menu item.

              CUSTOMIZATION OF GIDS:
              GIDS allows the user to specify some defaults in the .Xdefaults
              or .Xresources file. The parameters which can be set are the
              following:

              gids.calibration_font:     font for calibration window.
              gids.help_file:            path to help file.
              gids.help_font:            font for help window.
              gids.info_font:            font for info windows.
              gids.menu_font:            font for menu windows.
              gids.display_height:       default height of display window.
              gids.display_width:        default width of display window.
              gids.start_lut:            color lookup table at program start
                                         (use name of LUT as from color menu).
              gids.scaling:              color transfer funtion as program
                                         start (use name of scaling funtion
                                         as from scale menu).

              gids.inverse_scaling:      True for negative slope of color
                                         transfer function.
              gids.autoraise:            Raise window when new image or plot
                                         is arriving (default False).

(1)           ZOOM

              Pressing ZOOM will display record number zero, enter zoom mode
              and activate the zoom-menu.

              Mouse Button Control in display window:
              Button1     zoom in the image at the current position.
              Button2     re-centre the image at the current position.
              Button3     zoom out the image at the current position.

              Zooming can only affect the image which has been loaded into
              GIDS most recently. This image is always recorded in record 0.

(1.0)         ZOOM RESET

              Reset the image (i.e. no zoom factors).

(1.1)         ZOOM IN

              Zooms in the image by a factor 2.

(1.2)         ZOOM OUT

              Zooms out the image by a factor 2.

(1.3)         ZOOM IN-X

              Zooms in the image along the X axis by a factor 2.

(1.4)         ZOOM OUT-X

              Zooms out the image along the X axis by a factor 2.

(1.5)         ZOOM IN-Y

              Zooms in the image along the Y axis by a factor 2.

(1.6)         ZOOM OUT-Y

              Zooms out the image along the Y axis by a factor 2.

(1.7)         ZOOM RETURN

              Returns to main menu.

(2)           COLOR

              Pressing COLOR will enter color mode and activate the color-menu.

              Mouse Button Control in display window:
              Pressing any button and moving it in the display window will
              modify the color transfer function. Moving in the X direction
              will modify the slope of the transfer function (slope is zero at
              the right-most cursor position). Moving in the Y direction will
              modify the offset of the transfer function (largest offset at
              the upper-most cursor position).

              Pressing button3 will activate the monochrome LUT.
              Pressing button1 will change back to the original LUT active
              before pressing button3.
              Pressing button 2 will not change the LUT.

              The color-menu will allow the user to choose a color lookup table
              (LUT).

              NOTE: Color mode can always be activated by pressing the
              control key <Ctrl> on your keyboard and then the wanted
              mouse button.

(2.0)         COLOR RESET

              Resets the color transfer function to its initial shape (i.e.
              slope one and zero offset).

(2.1)         COLOR MONO

              Selects the monochrome LUT.

(2.2)         COLOR RAINBOW

              Selects the rainbow LUT.

(2.3)         COLOR RONEKERS

              Selects the LUT developped by Ron Ekers on the pdp 11/70.

(2.4)         COLOR STAIRCASE

              Selects the staircase LUT.

(2.5)         COLOR USER

              Selects the LUT which has been loaded by the user. A special
              program is needed for this (i.e. GIPSY program colour).

(2.6)         COLOR OTHER

              Gives another set of LUTs to choose from.

(2.7)         COLOR RETURN

              Returns to main menu, exits color mode.

(2.6.0)       COLOR OTHER BACKGROUND

              Selects the background LUT.

(2.6.1)       COLOR OTHER HEAT

              Selects the heat LUT.

(2.6.2)       COLOR OTHER ISOPHOT

              Selects the isophot LUT.

(2.6.3)       COLOR OTHER MOUSSE

              Selects the mousse LUT.

(2.6.4)       COLOR OTHER RANDOM

              Selects the random LUT, not at random.

(2.6.5)       COLOR OTHER RGB

              Selects the red-green-and-blue LUT.

(2.6.6)       COLOR OTHER SMOOTH

              Selects the smooth LUT.

(2.6.7)       COLOR OTHER RETURN

              Goes up one menu.

(3)           SCALE

              Pressing SCALE will enter color mode and activate the scale-menu.

              The scale-menu will allow the user to choose a particular color
              transfer function.

(3.0)         SCALE INVERSE

              Pressing INVERSE will make the slope of the transfer function
              negative.

(3.1)         SCALE LINEAR

              Pressing LINEAR will reset the color transfer function to its
              original shape (i.e. linear).

(3.2)         SCALE LOG

              Pressing LOG will make the color transfer function logarithmic
              (i.e. more color resolution at the lower end).

(3.3)         SCALE EXP

              Pressing EXP will make the color transfer function exponential
              (i.e. more color resolution at the upper end).

(3.4)         SCALE SQRT

              Pressing SQRT will make the color transfer function square root.

(3.5)         SCALE SQUARE

              Pressing SQRT will make the color transfer function square.

(3.6)         SCALE BLANK

              Entering 3 numbers between 0.0 and 1.0 will set resp. the
              Red, Green and Blue intensity for the BLANK data values.

(3.7)         SCALE RETURN

              Returns to main menu, exits color mode.

(4)           WEDGE

              Pressing WEDGE will enter wedge mode combined with color mode
              and activate the wedge menu.

              Mouse Button Control in wedge window:
              Pressing any button and moving it in the wedge window will
              change the color of the correspondig data values to the color
              selected with the wedge menu.

(4.0)         WEDGE RESET

              Pressing RESET will remove all wedge colors selected by the
              user.

(4.1)         WEDGE RED

              Pressing RED will select a red wedge color.

(4.2)         WEDGE GREEN

              Pressing GREEN will select a green wedge color.

(4.3)         WEDGE BLUE

              Pressing BLUE will select a blue wedge color.

(4.4)         WEDGE YELLOW

              Pressing YELLOW will select a yellow wedge color.

(4.5)         WEDGE PURPLE

              Pressing PURPLE will select a purple wedge color.

(4.6)         WEDGE WHITE

              Pressing WHITE will select a white wedge color.

(4.7)         WEDGE RETURN

              Returns to main menu, exits wedge mode.

(5)           CURSOR

              Pressing CURSOR will enable cursor mode. In cursor mode
              the user can send grid positions or the lower-left and
              upper-right grids of a box to an accepting program (usually
              hermes).

              Mouse Button Control in display window.
              Pressing button1 in the display window will send the
              corresponding grid position to hermes.
              Pressing button2 and moving the cursor allows the user to draw
              a box on the display window. Pressing button2 again will allow
              the user to modify this box. Pressing button3 inside the box will
              send the box grids to hermes and remove the box, pressing
              button3 outside the box will only remove it.

(6)           ETCETERA

              Pressing ETCETERA will activate the etcetera menu, where some
              extra options are made available, like manipulating recordings and
              graphics planes. Also the quit and resize button are programmed
              here.

(6.0)         ETCETERA QUIT

              Pressing QUIT will allow the user to quit via the quit menu.

(6.0.0)       ETCETERA QUIT NO

              Don't quit and go back to etcetera menu.

(6.0.1)       ETCETERA QUIT YES

              Gids is stopped graciously. Always exit via this menu item.

(6.1)         ETCETERA GRAPHICS

              Pressing GRAPHICS will open an extra window for plotting.
              Without the extra plot window plot output will appear in the
              display window.
              Pressing GRAPHICS again will remove the extra plot window
              and any plot output will again be redirected to the display
              window.

              The default plot window size is 512 x 512. The user can specify
              the size by typing the wanted X size and Y size (those numbers
              should appear in the lower info window) and then pressing
              GRAPHICS.

(6.2)         ETCETERA RECORDING

              Pressing RECORDING will activate the recording menu, where
              the user can manipulate the recorded images.

              The third info window displays information about the recordings
              (REC= MR, NR, CR), where MR is the maximum number of recordings
              the server can handle, NR the number of images recorded and
              CR the current record number in the display window. Record 0 is
              special since here always the last image loaded into GIDS is
              stored. It is only possible to (de)zoom this image.

(6.2.0)       ETCETERA RECORDING RECORD

              Typing the wanted record number and pressing RECORD will
              record the image currently in the display window.

(6.2.1)       ETCETERA RECORDING VIEW

              Pressing VIEW will enter view mode which will allow the user to
              manually loop through the recorded images.

              Mouse Button Control in display window:
              Pressing button1 in the display window will show the previous
              recording, pressing button3 will display the next recording.

              When pressing VIEW, the default action is to show the current
              record number. Any other recording can be viewed by typing its
              record number and then pressing VIEW.

(6.2.2)       ETCETERA RECORDING BLINK

              This option will allow the user to blink two images.

              Typing one or two record numbers and pressing BLINK will
              enter blink mode. If only one record number was typed, record 0
              will be used for the first recording and the record number
              entered for the second.

              Mouse Button Control in display window:
              Pressing any button will set the time interval between blinks.
              Blink frequency increases to the left.

(6.2.3)       ETCETERA RECORDING LOOP

              This option will allow the user to make a movie of
              the recorded images.

              Typing the wanted record numbers and pressing LOOP will
              enter loop mode.

              Mouse Button Control in display window:
              Button1: Loop backward.
              Button2: set the loop speed. Speed increases to the left.
              Button3: Loop forward.

              NOTE: Some programs already define the loop for the user
              (i.e. the GIPSY program VIEW with RECORD=Y). The user
              then only needs to press LOOP and does not have to specify
              the record numbers.

(6.2.4)       ETCETERA RECORDING SPLIT

              This option will allow the user to view two images in split
              screen.

              Typing one or two record numbers and pressing SPLIT will
              enter split screen mode. In split screen mode gids displays
              in different parts of the display window the two recordings
              entered by the user. If only one record number was typed,
              record zero will be used for the first recording and the
              record number entered for the second.

              Mouse Button Control in display window:
              Button1 will display to the left of the cursor the first
              recording entered by the user and to the right the second.
              Button2 will display above the cursor the second recording
              entered by the user and below the cursor the first.
              Button3 will display in the right-upper part above the cursor
              the second recording the user entered.

(6.2.5)       ETCETERA RECORDING REMOVE

              Typing a record number and pressing REMOVE will remove the
              recording.

(6.2.6)       ETCETERA RECORDING RETURN

              Returns to etcetera menu.

(6.3)         ETCETERA REGION

              Pressing REGION will display record number zero and activate
              the region menu.

              This option is only of use in combination with a special
              program which can do something useful with the graphics data
              (i.e. GIPSY program BLOT for masking data).

              The fourth info window displays information about the graphics
              planes (MASK=   M ON=  N), where M is the current mask of
              graphics planes on which the selected operation (i.e. DEFINE)
              is valid and N the mask for the planes which are displayed.
              GIDS has four graphics planes which are numbered 1, 2, 4 and 8.
              If the user wants more than one graphics plane simply add the
              numbers.

              Graphics operations only work on the last image read into GIDS.

(6.3.0)       ETCETERA REGION DEFINE

              Pressing DEFINE will allow the user to start defining a region
              on top of record zero. The graphics planes which will be modified
              are shown in the fourth info window (MASK=).

              Mouse Button Control in display window:
              Pressing Button1 for the first time will start the region
              (a polygon) on that position. The following clicks on Button1
              will draw the outer boundaries of the polygon.
              Pressing Button2 will remove the last drawn boundaries.
              Pressing Button3 will close and fill the polygon.

(6.3.1)       ETCETERA REGION UNDEFINE

              Pressing UNDEFINE will alow the user to start defining a region
              which is to be removed from a previous defined region. The
              graphics planes which will be modified are shown in the fourth
              info window (MASK=).

              Mouse Button Control in display window:
              Pressing Button1 for the first time will start the region
              (a polygon) to be removed on that position. The following clicks
              on Button1 will draw the outer boundaries of the polygon to be
              removed.
              Pressing Button2 will remove the last drawn boundaries.
              Pressing Button3 will close and remove the polygon.

(6.3.2)       ETCETERA REGION ON

              Typing a new ON mask (displayed in the fourth info window)
              will display the graphics data in the planes which match
              the new ON mask.

(6.3.3)       ETCETERA REGION CLEAR

              Pressing CLEAR will clear the graphics planes which match
              the MASK shown in the fourth info window.

              The user can modify this mask by typing the wanted mask and
              pressing CLEAR.

(6.3.4)       ETCETERA REGION READY

              Press READY when finished with (un)defining regions. This
              informs the program which is interested in the graphics plane
              data to read it.

(6.3.5)       ETCETERA REGION RETURN

              Returns to etcetera menu.

(6.4)         ETCETERA RESIZE

              Typing the wanted width and height in pixels of the display
              window and pressing RESIZE will resize the display window.

              Note that all previous recordings will be lost!

(6.5)         ETCETERA SLICE

              Pressing SLICE will enter slice mode and activate the
              slice menu.

              The slice menu allows the user to display cuts along a principal
              axis through all recordings. The xslice window (created with
              the X-SLICE option) displays an X-cut through all selected
              recordings. The yslice window (created with the Y-SLICE option)
              displays a Y-cut through all selected recordings.

              Mouse Button Control in display window:
              Pressing any button will show in the xslice window (the window
              below the display window) an X-cut through all selected
              recordings at the present x position, and in the yslice window
              (the window to the right of the display window) a Y-cut through
              all selected recordings at the present y position.

              Mouse Button Control in xslice window:
              Pressing any button will show in the display window the
              corresponding recorded image and in the yslice window the
              corresponding Y-cut.

              Mouse Button Control in yslice window:
              Pressing any button will show in the display window the
              corresponding recorded image and in the xslice window the
              corresponding X-cut.

              Note: Pressing button2 in any of the above windows will
              display a cursor to indicate the cross-sections. Any other
              button will remove the cursor.

(6.5.0)       ETCETERA SLICE X-SLICE

              Entering the wanted recorded images and then pressing X-SLICE
              will create the xslice window below the display window. The
              xslice window displays an X-cut through the selected recordings
              at the position of the cursor in the display window when
              pressing a mouse button.

              Pressing X-SLICE with no recordings selected will remove the
              xslice window.

(6.5.1)       ETCETERA SLICE Y-SLICE

              Entering the wanted recorded images and then pressing Y-SLICE
              will create the yslice window to the right of the display
              window. The yslice window displays a Y-cut through the selected
              recordings at the position of the cursor in the display window
              when pressing a mouse button.

              Pressing Y-SLICE with no recordings selected will remove the
              yslice window.

(6.5.2)       ETCETERA SLICE ZOOM-X

              Pressing ZOOM-X will double every X-cut in the xslice window.

(6.5.3)       ETCETERA SLICE ZOOM-Y

              Pressing ZOOM-Y will double every Y-cut in the yslice window.

(6.5.4)       ETCETERA SLICE DEZOOM-X

              Pressing DEZOOM-X will remove every second X-cut in the xslice
              window.

(6.5.5)       ETCETERA SLICE DEZOOM-Y

              Pressing DEZOOM-Y will remove every second Y-cut in the yslice
              window.

(6.5.6)       ETCETERA SLICE RETURN

              Returns to etcetera menu and leaves slice mode.

(6.6)         ETCETERA RETURN

              Returns to main menu.
.fi
.SH AUTHOR
K.G. Begeman
.SH SEE ALSO
fitsgids(1NEMO), mirgids(1NEMO)
.SH HISTORY
.nf
Jul 17, 1991: KGB, Document created.
Aug 18, 1993: KGB, Document improved.
Jul 22, 1994: PJT, Unix manual page
.fi