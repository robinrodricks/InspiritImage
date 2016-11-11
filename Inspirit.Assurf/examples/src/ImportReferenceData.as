package
{
	import ru.inspirit.surf.ASSURF;
	import ru.inspirit.surf.ReferenceInfo;
	import ru.inspirit.surf.ar.ARCamera3D;
	import ru.inspirit.surf.arsupport.demo.away3dlite.World3D;
	import ru.inspirit.surf.image.AutoImageProcessor;
	import ru.inspirit.surf.image.QuasimondoImageProcessor;

	import flash.display.Bitmap;
	import flash.display.BitmapData;
	import flash.display.Sprite;
	import flash.display.StageScaleMode;
	import flash.events.Event;
	import flash.events.TimerEvent;
	import flash.geom.Rectangle;
	import flash.media.Camera;
	import flash.media.Video;
	import flash.text.TextField;
	import flash.utils.ByteArray;
	import flash.utils.Endian;
	import flash.utils.Timer;
	import flash.utils.getTimer;

	/**
	 * @author Eugene Zatepyakin
	 */
	 
	[SWF(width='640',height='590',frameRate='25',backgroundColor='0xFFFFFF')]
	
	public class ImportReferenceData extends Sprite
	{
		[Embed('../assets/points_low.ass', mimeType='application/octet-stream')] protected static const ref_data_ass:Class;
		
		protected var myview:Sprite;
		protected var _txt:TextField;
		protected var camBmp:Bitmap;
		
		protected var _cam:Camera;
		protected var _video:Video;
		protected var _cambuff:BitmapData;
		protected var _cambuff_rect:Rectangle;
		
		public const surf:ASSURF = new ASSURF();
		public const imgQProc:QuasimondoImageProcessor = new QuasimondoImageProcessor(new Rectangle(), false);
		public const imgAProc:AutoImageProcessor = new AutoImageProcessor(new Rectangle());
		public const surfTimer:Timer = new Timer(1000/20);
		
		public var arCamera:ARCamera3D = new ARCamera3D();
		public var world3d:World3D;
		
		public var refID_0:int;
		
		public function ImportReferenceData()
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
			
			// ASSURF setup
			// first method you should call
			surf.init(ASSURF.DETECT_PRECISION_LOW, 350, 10000, 1);
			
			// make ASSURF detect region of interest automatically
			//surf.autoDetectROI = true;
			
			// add ByteArray object as reference
			var ref_data:ByteArray = ByteArray(new ref_data_ass());
			ref_data.endian = Endian.LITTLE_ENDIAN;
			ref_data.uncompress();
			ref_data.position = 0;
			surf.importReferenceData(ref_data);
			refID_0 = 0;
			
			// finalize setup by analyzing all references
			// you cant add anything after it
			surf.buildRefIndex();
			
			// here is 2 different image pre-processors you can use to refine camera image
			imgAProc.imageRect = imgQProc.imageRect = _cambuff.rect;
			surf.imageProcessor = imgAProc; // smoother
			//surf.imageProcessor = imgQProc; // harder
			
			// setup input image dimensions
			surf.setup(640, 480);
			
			addChild(myview);
			
			// init Away3DLite objects
			init3d();
			
			surfTimer.addEventListener(TimerEvent.TIMER, process);
			addEventListener(Event.ENTER_FRAME, render3d);
			surfTimer.start();
		}
		
		protected function process(e:Event = null):void
		{
			var t:int = getTimer();
			var ref:ReferenceInfo;

			ref = surf.detectSingleObject(_cambuff, refID_0, ASSURF.GET_3DPOSE);
			
			if(ref.matchedPointsCount > 4 && ref.transformError < 5)
			{
				world3d.airplane.setTransform(ref.transform.clone(), ref.transformError);
			} 
			else 
			{
				world3d.airplane.lost();
			}
			
			_txt.text = surf.debug();
			_txt.appendText('\npose err: ' + ref.transformError);
			_txt.appendText('\ndone in: ' + String(getTimer()-t) + 'ms');
		}
		
		protected function render3d(e:Event = null):void
		{
			_cambuff.draw(_video);
			
			world3d.render();
		}

		protected function initCamera(w:int = 640, h:int = 480, fps:int = 25):void
		{
			_cambuff = new BitmapData( w, h, false, 0x0 );
			_cam = Camera.getCamera();
			_cam.setMode( w, h, fps, true );
			
			_video = new Video( _cam.width, _cam.height );
			_video.attachCamera( _cam );
		}
		
		protected function init3d():void
		{
            world3d = new World3D(surf, arCamera);            
			world3d.initPlane();
			
			addChild(world3d);
		}
	}
}
