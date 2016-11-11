package demo.pyrFlowLK
{
	import apparat.math.FastMath;
	import apparat.memory.Memory;

	import ru.inspirit.image.ImagePyramid;
	import ru.inspirit.image.flow.OpticalFlowPyrLK;

	import flash.display.Bitmap;
	import flash.display.BitmapData;
	import flash.display.Sprite;
	import flash.display.StageScaleMode;
	import flash.events.Event;
	import flash.events.MouseEvent;
	import flash.geom.Matrix;
	import flash.geom.Rectangle;
	import flash.media.Camera;
	import flash.media.Video;
	import flash.text.TextField;
	import flash.utils.ByteArray;
	import flash.utils.Endian;
	import flash.utils.getTimer;

	/**
	 * @author Eugene Zatepyakin
	 */
	[SWF(width='640',height='590',frameRate='25',backgroundColor='0xFFFFFF')]
	
	public final class DemoPyrFlowLK extends Sprite
	{		
		protected var myview:Sprite;
        public static var _txt:TextField;
        protected var camBmp:Bitmap;
        
        protected var _cam:Camera;
        protected var _video:Video;
        protected var _cambuff:BitmapData;
        protected var _cambuff_rect:Rectangle;
        protected var _cam_mtx:Matrix;
        
        public const ram:ByteArray = new ByteArray();
        public const flowPyrLK:OpticalFlowPyrLK = new OpticalFlowPyrLK();
        public const imgPyr:ImagePyramid = new ImagePyramid(3);
        public const prevPyr:ImagePyramid = new ImagePyramid(3);
        
        public var trackPoints:Vector.<TrackPoint> = new Vector.<TrackPoint>();
        
        public var LK_WIN:int = 21;
        
        public var frame:int = 0;
        public var frame_prev:int = 0;
        public var frame_ms:int = 0;
		
		public function DemoPyrFlowLK()
		{
			if(stage) init();
			else addEventListener(Event.ADDED_TO_STAGE, init);
		}
		
		protected function init(e:Event = null):void
		{
			removeEventListener(Event.ADDED_TO_STAGE, init);
			//
			stage.scaleMode = StageScaleMode.NO_SCALE;

            myview = new Sprite();
            
            // debug test field
            _txt = new TextField();
            _txt.autoSize = 'left';
            _txt.width = 300;
            _txt.x = 5;
            _txt.y = 480;                   
            myview.addChild(_txt);
            
            // web camera initiation
            initCamera(640, 480, 25);
            camBmp = new Bitmap(_cambuff);                  
            myview.addChild(camBmp);
            
			var chunkLK:int = flowPyrLK.calcRequiredChunkSize(LK_WIN, 3);
			var chunkPyr:int = imgPyr.calcRequiredChunkSize(640, 480);
			var chunkPyr2:int = prevPyr.calcRequiredChunkSize(640, 480);

			ram.endian = Endian.LITTLE_ENDIAN;
			ram.length = chunkLK + chunkPyr + chunkPyr2 + 1024;
			ram.position = 0;
			Memory.select(ram);

			flowPyrLK.setup( 1024, 640, 480, LK_WIN, 3 );
			
			imgPyr.setup( 640, 480, 1024 + chunkLK );
			prevPyr.setup( 640, 480, 1024 + chunkLK + chunkPyr );
			imgPyr.srcImage = prevPyr.srcImage = _cambuff;
            
            addChild(myview);
            
			addEventListener(Event.ENTER_FRAME, render);
			stage.addEventListener(MouseEvent.CLICK, addPoint);
		}
		
		protected function render(e:Event = null):void
		{
			var t:int = getTimer();
			
			_cambuff.draw(_video, _cam_mtx);
			
			imgPyr.updateMem( true, 3, 3 );
			//imgPyr.mem[0].render( _cambuff );
			//imgPyr.mem[1].render( _cambuff );
			//imgPyr.mem[2].render( _cambuff );
			
			var pyrt:int = getTimer() - t;
			
			var n:int = trackPoints.length;
			var newPoints:Vector.<Number> = new Vector.<Number>(n<<1);
			var prevPoints:Vector.<Number> = new Vector.<Number>(n<<1);
			var status:Vector.<int> = new Vector.<int>(n);
			var pp:TrackPoint;
			var fx:Number, fy:Number;
			for(var i:int = 0; i < n; ++i)
			{
				pp = trackPoints[i];
				newPoints[i<<1] = pp.x + pp.vx;
				newPoints[((i<<1)+1)|0] = pp.y + pp.vy;
				prevPoints[i<<1] = pp.x;
				prevPoints[((i<<1)+1)|0] = pp.y;

				pp.tracked = false;
			}
			
			t = getTimer();
			
			flowPyrLK.calcOpticalFlow(
										Vector.<int>([prevPyr.mem[0].ptr, prevPyr.mem[1].ptr, prevPyr.mem[2].ptr]), 
										Vector.<int>([imgPyr.mem[0].ptr, imgPyr.mem[1].ptr, imgPyr.mem[2].ptr]), 
										prevPoints, newPoints, n, status, null, OpticalFlowPyrLK.INITIAL_GUESSES, 20, 0.025);
										
			var flowt:int = getTimer() - t;
							
			for(i = 0; i < n; ++i)
			{
				if(status[i])
				{
					fx = newPoints[i<<1];
					fy = newPoints[((i<<1)+1)|0];
					pp = trackPoints[i];
					pp.vx = (fx - pp.x) * 0.75;
					pp.vy = (fy - pp.y) * 0.75;
					pp.x = fx;
					pp.y = fy;
					pp.tracked = true;
				}
			}

			plotPoints(trackPoints, n);
			
			// swap images
			imgPyr.swap(prevPyr);
			
			_txt.text = '3 Level Image Pyramid | 21px Patch Size | 20 Iterations per point | '+ n +' points\n';
			if(getTimer() - frame_ms >= 1000)
			{
				_txt.appendText( frame + '/25 fps\nimage pyramid: ' + pyrt + 'ms\nopt flow lk: ' + flowt + 'ms' );
				frame_prev = frame;
				frame_ms = getTimer();
				frame = 0;
			}
			else
			{
				frame++;
				
				_txt.appendText( frame_prev + '/25 fps\nimage pyramid: ' + pyrt + 'ms\nopt flow lk: ' + flowt + 'ms' );
			}
		}
		
		public function addPoint(e:Event = null):void
		{
			var n:int = trackPoints.length;
			var p:TrackPoint;
			var nx:Number = mouseX;
			var ny:Number = mouseY;
			for(var i:int = 0; i < n; ++i)
			{
				p = trackPoints[i];
				var dx:Number = p.x - nx;
				var dy:Number = p.y - ny;
				if(dx*dx + dy*dy < 100)
				{
					trackPoints.splice(i, 1);
					return;
				}
			}
			
			p = new TrackPoint();
			p.x = nx;
			p.y = ny;
			p.vx = p.vy = 0;
			trackPoints.push(p);
		}
		
		protected function plotPoints(pts:Vector.<TrackPoint>, n:int):void
		{
			var px:int, py:int;
			var col:uint;
			
			_cambuff.lock();
			
			for(var i:int = 0; i < n; ++i)
			{
				px = FastMath.rint(pts[i].x);
				py = FastMath.rint(pts[i].y);
				col = pts[i].tracked ? 0x00FF00 : 0xFF0000;
				
				_cambuff.setPixel(px, py, col);
				_cambuff.setPixel(px+1, py, col);
				_cambuff.setPixel(px-1, py, col);
				_cambuff.setPixel(px, py+1, col);
				_cambuff.setPixel(px, py-1, col);
			}
			_cambuff.unlock( _cambuff_rect );
		}
		
		protected function initCamera(w:int = 640, h:int = 480, fps:int = 25):void
        {
            _cambuff = new BitmapData( w, h, false, 0x0 );
            _cam = Camera.getCamera();
            _cam.setMode( w, h, fps, true );

			_cambuff_rect = _cambuff.rect;
			_cam_mtx = new Matrix(-1, 0, 0, 1, w);
            
            _video = new Video( _cam.width, _cam.height );
            _video.attachCamera( _cam );
        }
	}
}
