package 
{
	import ru.inspirit.surf.ASSURF;
	import ru.inspirit.surf.ReferenceInfo;
	import ru.inspirit.surf.image.AutoImageProcessor;
	import ru.inspirit.surf.image.QuasimondoImageProcessor;

	import flash.display.Bitmap;
	import flash.display.BitmapData;
	import flash.display.Graphics;
	import flash.display.Shape;
	import flash.display.Sprite;
	import flash.display.StageScaleMode;
	import flash.events.Event;
	import flash.events.TimerEvent;
	import flash.geom.Point;
	import flash.geom.Rectangle;
	import flash.media.Camera;
	import flash.media.Video;
	import flash.text.TextField;
	import flash.utils.Timer;
	import flash.utils.getTimer;
	
	/**
	 * @author Eugene Zatepyakin
	 */

	[SWF(width='640',height='590',frameRate='25',backgroundColor='0xFFFFFF')]

	public final class Homography extends Sprite
	{
		[Embed(source = '../assets/graffiti_800.jpg')] private static const ref_a:Class;
		
		protected const iref : BitmapData = Bitmap(new ref_a()).bitmapData;
		
		protected var myview:Sprite;
		protected var _txt:TextField;
		protected var camBmp:Bitmap;
		protected var outline:Shape;
		
		protected var _cam:Camera;
		protected var _video:Video;
		protected var _cambuff:BitmapData;
		protected var _cambuff_rect:Rectangle;
		
		public const surf:ASSURF = new ASSURF();
		public const imgQProc:QuasimondoImageProcessor = new QuasimondoImageProcessor(new Rectangle(), false);
		public const imgAProc:AutoImageProcessor = new AutoImageProcessor(new Rectangle());
		public const surfTimer:Timer = new Timer(1000/20);
		
		public var refID_0:int;

		public function Homography()
		{			
			if(stage) init();
			else addEventListener(Event.ADDED_TO_STAGE, init);
		}

		protected function init():void
		{
			removeEventListener(Event.ADDED_TO_STAGE, init);
			
			stage.scaleMode = StageScaleMode.NO_SCALE;

			myview = new Sprite();
			
			// debug test field
			_txt = new TextField();
			_txt.autoSize = 'left';
			_txt.width = 300;
			_txt.x = 10;
			_txt.y = 480;			
			myview.addChild(_txt);
			
			// web camera initiation
			initCamera(640, 480, 25);
			camBmp = new Bitmap(_cambuff);			
			myview.addChild(camBmp);
			
			outline = new Shape();
			myview.addChild(outline);
			
			// ASSURF setup
			// first method you should call
			surf.init(ASSURF.DETECT_PRECISION_LOW, 300, 10000, 1);
			
			// make ASSURF detect region of interest automatically
			surf.autoDetectROI = true;
			
			// add bitmapData object as reference
			refID_0 = surf.addRefObject(iref, 5, 1800, false);
			
			// finalize setup by analyzing all references
			// you cant add anything after it
			surf.buildRefIndex();
			
			// here is 2 different image pre-processors you can use to refine camera image
			imgAProc.imageRect = imgQProc.imageRect = _cambuff.rect;
			surf.imageProcessor = imgAProc;
			//surf.imageProcessor = imgQProc;
			
			// setup input image dimensions
			surf.setup(640, 480);
			
			addChild(myview);
			
			surfTimer.addEventListener(TimerEvent.TIMER, process);
			addEventListener(Event.ENTER_FRAME, drawCamera);
			surfTimer.start();
		}
		
		protected function process(e:Event = null):void
		{
			var t:int = getTimer();
			var ref:ReferenceInfo;
			var gfx:Graphics = outline.graphics;
			
			ref = surf.detectSingleObject(_cambuff, refID_0, ASSURF.GET_HOMOGRAPHY);
			
			gfx.clear();
			gfx.lineStyle(2, 0x00FF00);
			
			if(ref.matchedPointsCount > 4) 
			{
				var pt0 : Point = ref.homography.projectPoint(new Point(0, 0));
				var pt1 : Point = ref.homography.projectPoint(new Point(iref.width, 0));
				var pt2 : Point = ref.homography.projectPoint( new Point(iref.width, iref.height) );
				var pt3 : Point = ref.homography.projectPoint( new Point(0, iref.height) );
				gfx.moveTo(pt0.x, pt0.y);
				gfx.lineTo(pt1.x, pt1.y);
				gfx.lineTo(pt2.x, pt2.y);
				gfx.lineTo(pt3.x, pt3.y);
				gfx.lineTo(pt0.x, pt0.y);
			}
			
			_txt.text = surf.debug();
			_txt.appendText('\ndone in: ' + String(getTimer()-t) + 'ms');
		}
		
		protected function drawCamera(e:Event):void
		{
			_cambuff.draw(_video);
		}

		protected function initCamera(w:int = 640, h:int = 480, fps:int = 25):void
		{
			_cambuff = new BitmapData( w, h, false, 0x0 );
			_cam = Camera.getCamera();
			_cam.setMode( w, h, fps, true );
			
			_video = new Video( _cam.width, _cam.height );
			_video.attachCamera( _cam );
		}
	}
}
