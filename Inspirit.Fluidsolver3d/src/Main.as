package
{
	import ru.inspirit.fluid.FluidRenderer;
	import ru.inspirit.fluid.FluidSolver3D;
	import ru.inspirit.fluid.Particle3DManager;
	import ru.inspirit.fluid.emitter.IFluidEmitter;
	import ru.inspirit.fluid.emitter.FluidForceEmitter;
	import ru.inspirit.fluid.emitter.FluidColorEmitter;
	import ru.inspirit.fluid.emitter.FluidMovingEmitter;

	import com.bit101.components.Label;

	import org.audiofx.mp3.MP3FileReferenceLoader;
	import org.audiofx.mp3.MP3SoundEvent;

	import flash.display.Bitmap;
	import flash.display.BitmapData;
	import flash.display.BlendMode;
	import flash.display.Sprite;
	import flash.display.StageDisplayState;
	import flash.display.StageScaleMode;
	import flash.events.ContextMenuEvent;
	import flash.events.Event;
	import flash.events.KeyboardEvent;
	import flash.events.MouseEvent;
	import flash.filters.BlurFilter;
	import flash.geom.ColorTransform;
	import flash.geom.Matrix;
	import flash.geom.Matrix3D;
	import flash.geom.PerspectiveProjection;
	import flash.geom.Point;
	import flash.geom.Rectangle;
	import flash.geom.Vector3D;
	import flash.media.Sound;
	import flash.media.SoundChannel;
	import flash.media.SoundMixer;
	import flash.net.FileFilter;
	import flash.net.FileReference;
	import flash.system.System;
	import flash.ui.ContextMenu;
	import flash.ui.ContextMenuItem;
	import flash.ui.Keyboard;
	import flash.utils.ByteArray;
	import flash.utils.getTimer;

	/**
	 * Fluid Solver test application
     * @author Eugene Zatepyakin
     */

	[SWF(width='780',height='460',frameRate='60',backgroundColor='0x000000')]

	public class Main extends Sprite
	{

		protected const origin:Point = new Point();
		protected const identity:Matrix = new Matrix();
		protected const blur:BlurFilter = new BlurFilter( 8, 8, 2 );
		protected const fade2alpha:ColorTransform = new ColorTransform( 1, 1, 1, .9);
		
		public const W:int = 780;
		public const H:int = 460;

		public const scrW:int = 460;
		public const scrH:int = 460;
		public const scrD:int = 460;

		public const solverW:int = 25;
		public const solverH:int = int( solverW * scrH / scrW );
		public const solverD:int = int( solverW * scrD / scrW );
		public const solver:FluidSolver3D = new FluidSolver3D(solverW, solverH, solverD);
		public const particles:Particle3DManager = new Particle3DManager(solver, 0, 0, scrW, scrH, 0, scrD);
		public const fluids:FluidRenderer = new FluidRenderer(solver);
		
		public const emitters:Vector.<IFluidEmitter> = new Vector.<IFluidEmitter>();

		public const proj:PerspectiveProjection = new PerspectiveProjection();
		public const tMatrix:Matrix3D = new Matrix3D();
		public var pMatrix:Matrix3D;
		public var focalLength:Number;
		public var focalLength2:Number;
		
		public const fluids_bmp:BitmapData = new BitmapData(W*0.25, H*0.25, false, 0x00);
		public const fluids_bmp2:BitmapData = new BitmapData(W, H, false, 0x00);
		public const particles_bmp:BitmapData = new BitmapData(W, H, true, 0);
		public const particles_canvas:Bitmap = new Bitmap(particles_bmp, 'never', false);
		public const fluids_canvas:Bitmap = new Bitmap(fluids_bmp, 'always', true);
		public const fluids_canvas2:Bitmap = new Bitmap(fluids_bmp2, 'never', false);
		
		protected const sparkleScale:int = 2; 
		protected const sparkle_bmp:BitmapData = new BitmapData(scrW/sparkleScale, scrH/sparkleScale, true, 0x0);
		protected const sparkle_canvas:Bitmap = new Bitmap(sparkle_bmp, 'never', true);
		protected const sparkleDrawMatrix:Matrix = new Matrix(1/sparkleScale, 0, 0, 1/sparkleScale, 0, 0);
		
		public var isLive:Boolean = false;
		public var sound:Sound;
		public var sndCh:SoundChannel;
		
		protected var spectrum:ByteArray = new ByteArray();
		protected var spinDir:Number = 1;
		protected var angle:Number = 0;
		protected var fx:Number = 0;
		protected var fy:Number = 0;
		protected var fz:Number = 0;
		protected var dx:Number = 1;
		protected var dy:Number = -1;
		protected var dz:Number = 1;

		protected var _dVel:Number = 0;
		protected var _dPhase:Number = 0;
		protected var solve_t:int = 0;
		protected var render_t:int = 0;
		protected var count_n:int = 0;
		protected var _timer:uint;
		protected var _fps:uint;
		protected var _ms:uint;
		protected var _ms_prev:uint;
		protected var frameCount:uint = 0;
		protected var info:Label;

		public function Main()
		{
			if (stage) init();
			else addEventListener(Event.ADDED_TO_STAGE, init);
		}

		protected function init(e:Event = null):void
		{
			removeEventListener(Event.ADDED_TO_STAGE, init);
			initStage();
			
			proj.fieldOfView = 75;
			focalLength = proj.focalLength;
			pMatrix = proj.toMatrix3D();
			
			var falf_tan:Number = Math.tan(75 / 2 * Math.PI / 180);
			focalLength = particles_bmp.width/2/falf_tan;
			focalLength2 = fluids_bmp.width/2/falf_tan;
			
			//fluids_canvas.width = scrW;
			//fluids_canvas.height = scrH;
			
			fluids_canvas.scaleX = fluids_canvas.scaleY = 4;
			fluids_canvas.x = (W - fluids_canvas.width) >> 1;
			fluids_canvas.y = (H - fluids_canvas.height) >> 1;
			
			//fluids_canvas2.scaleX = fluids_canvas2.scaleY = 8;
			//fluids_canvas2.x = (scrW - fluids_canvas2.width) >> 1;
			//fluids_canvas2.y = (scrH - fluids_canvas2.height) >> 1;
			
			sparkle_canvas.scaleX = sparkle_canvas.scaleY = sparkleScale;
			sparkle_canvas.blendMode = BlendMode.ADD;
			particles_canvas.blendMode = BlendMode.ADD;
			
			addChild(fluids_canvas);
			//addChild(fluids_canvas2);
			
			addChild(particles_canvas);
			addChild(sparkle_canvas);
			
			info = new Label(this, 10, 5);
			
			createEmitters();

			addEventListener(Event.ENTER_FRAME, render);
		}

		protected function render(e:Event):void
		{	       
			if(isLive)
			{
				computeSpectrum();     
			}
			
			var t:int = getTimer();
			
			solver.update();
			
			solve_t += getTimer() - t;
			t = getTimer();
			
			tMatrix.identity();
			
			if(isLive)
			{
				tMatrix.appendRotation( Math.cos(t * 0.00077) * 60, Vector3D.X_AXIS );	            
				tMatrix.appendRotation( 45+Math.cos(t * 0.00091) * 120, Vector3D.Y_AXIS );
			} else {
				tMatrix.appendRotation( (stage.mouseX - W * 0.5), Vector3D.Y_AXIS );
			}
			
			tMatrix.appendTranslation(0, 0, scrD*0.25);
			
			particles_bmp.lock();
			fluids_bmp.lock();
			sparkle_bmp.lock();
			
			particles_bmp.fillRect(particles_bmp.rect, 0x0);
			//fluids_bmp2.fillRect(fluids_bmp2.rect, 0x0);
			//fluids_bmp.fillRect(fluids_bmp.rect, 0x0);
			
			fluids.render(fluids_bmp, tMatrix, focalLength2);
			
			tMatrix.identity();				
			
			if(isLive)
			{
				tMatrix.appendRotation( Math.cos(t * 0.00077) * 60, Vector3D.X_AXIS );	            
				tMatrix.appendRotation( 45+Math.cos(t * 0.00091) * 120, Vector3D.Y_AXIS );
			} else 
			{
				tMatrix.appendRotation( (stage.mouseX - W * 0.5), Vector3D.Y_AXIS );
			}
			tMatrix.appendTranslation(0, 0, scrD);
			
			particles.render(particles_bmp, tMatrix, focalLength);
			
			fluids_bmp.applyFilter( fluids_bmp, fluids_bmp.rect, origin, blur );
			//fluids_bmp2.applyFilter( fluids_bmp2, fluids_bmp2.rect, origin, blur );

			if(sparkle_canvas.visible)
			{
	            sparkle_bmp.fillRect(sparkle_bmp.rect, 0x0);
	            sparkle_bmp.draw(particles_bmp, sparkleDrawMatrix);
			}
            
            sparkle_bmp.unlock(sparkle_bmp.rect);
            particles_bmp.unlock(particles_bmp.rect);
			fluids_bmp.unlock(fluids_bmp.rect);
			
			render_t += getTimer() - t;
			count_n++;
			
			var i:int;
			if(isLive)
			{
				for(i = 1; i < emitters.length; ++i) 
				{
					emitters[i].updateWithForce((fx+230)/scrW, (fy+230)/scrH, (fz+230)/scrD);
				}
			}
			else
			{
				for(i = 0; i < 4; ++i) 
				{
					emitters[i].update();
				}
			}
			
			countTime();
		}
		
		protected function createEmitters():void
		{
			var ind:int = -1;
			
			emitters[++ind] = new FluidForceEmitter(solver, particles, solverW >> 1, solverH >> 1, solverD >> 1);
			//emitters[++ind] = new FluidColorEmitter(solver, solverW >> 1, solverH >> 1, solverD >> 1);
			emitters[++ind] = new FluidMovingEmitter(solver, particles, solverW >> 1, solverH >> 1, solverD >> 1);
			emitters[++ind] = new FluidMovingEmitter(solver, particles, solverW >> 1, solverH >> 1, solverD >> 1);
			emitters[++ind] = new FluidMovingEmitter(solver, particles, solverW >> 1, solverH >> 1, solverD >> 1);
			
			emitters[++ind] = new FluidMovingEmitter(solver, particles, solverW >> 1, solverH >> 1, solverD >> 1);
			emitters[++ind] = new FluidMovingEmitter(solver, particles, solverW >> 1, solverH >> 1, solverD >> 1);
			emitters[++ind] = new FluidMovingEmitter(solver, particles, solverW >> 1, solverH >> 1, solverD >> 1);
		}
		
		protected function computeSpectrum():void 
		{
			SoundMixer.computeSpectrum(spectrum, true);
			
            var i:int = -1;
            var val:Number = 0;
            var forceMultiplier:Number = 9.5;
            var valx:Number = 0;
            var valy:Number = 0;
            var tembr:Number = 0;
            
            spectrum.position = 0;
            
            while (++i < 512)
            {
                
                val = spectrum.readFloat();
                if (i < 86)
                {
                    valx += val;
                }
                else if (i < 172)
                {
                    valy += val;
                }
                else
                {
                    tembr += val;
                }
            }
            
            var med:Number = int(valx + valy + tembr) / 3;
            
            if (med < 13 && Math.random() < 0.2)
            {
                spinDir *= -1;
            }
            
            solver.fadeSpeed = 1 - Math.min(1.001, med / 13);
            solver.solverIter = Math.max(2, Math.min(5, med));
            
            fx = Math.sin(angle) * valx * forceMultiplier;
            fy = (-Math.cos(angle)) * valy * forceMultiplier;
            fz = Math.sin(angle) * Math.cos(angle) * ((valx + valy) * 0.5) * forceMultiplier;
            
            //trace(fx, fy, fz);
            
            angle += valy / 100 * spinDir;
        }
        
        protected function toAutoMode(e:Event = null):void
        {
        	ContextMenuItem(contextMenu.customItems[1]).enabled = false;
        	isLive = false;
        	solver.fadeSpeed = 0.007;
            solver.solverIter = 4;
        	stopLive();
        }
        
        protected var loader:MP3FileReferenceLoader = new MP3FileReferenceLoader();
		protected var fileReference:FileReference = new FileReference();
		protected function loadMP3(e:Event = null):void
		{
			loader.removeEventListener(MP3SoundEvent.COMPLETE, mp3LoaderCompleteHandler);
			loader.addEventListener(MP3SoundEvent.COMPLETE, mp3LoaderCompleteHandler);
			
			fileReference.removeEventListener(Event.SELECT, fileReferenceSelectHandler);
			fileReference.addEventListener(Event.SELECT, fileReferenceSelectHandler);
			
			fileReference.browse([new FileFilter("mp3 files","*.mp3")]);
		}
		protected function fileReferenceSelectHandler(ev:Event):void
		{
			var name:String = fileReference.name;
			var ext	:String = name.substr(name.lastIndexOf("."), name.length);
			ext = ext.toLowerCase();
			
			if(ext == '.mp3') 
			{
				loader.getSound(fileReference);
			}
		}
		protected function mp3LoaderCompleteHandler(ev:MP3SoundEvent):void
		{
			if(isLive) stopLive();
			
			ContextMenuItem(contextMenu.customItems[1]).enabled = true;
			
			for(var i:int = 0; i < emitters.length; ++i) 
			{
				emitters[i].setRandomPosition();
			}
			
			sound = ev.sound;
			sndCh = sound.play();
			
			if(!isLive) 
			{
				isLive = true;
			}
		}
		
		protected function stopLive(e:Event = null):void
		{
			sndCh.stop();
			try {
				sound.close();
			}catch (e:Error) {}
			
			sndCh = null;
			sound = null;
		}
		
		protected function countTime():void
		{
			_timer = getTimer();
			
			if( _timer - 1000 >= _ms_prev )
			{
				_ms_prev = _timer;
				
				info.text = 'FPS: ' + _fps + '\tMEM: ' + Number((System.totalMemory * 0.000000954).toFixed(3))
					+'\tSOLVE TIME: ' + int(solve_t / count_n + 0.5) + 'ms\tRENDER TIME: ' + int(render_t / count_n + 0.5) + 'ms'
					+'\nSOLVER SIZE: ' + solverW + 'x' + solverH + 'x' + solverD + '\tPARTICLES: ' + particles.numUsedParticles;
				count_n = solve_t = render_t = _fps = 0;
			}
			
			_fps ++;
			_ms = _timer;
		}
		
		protected function toggleSparkles(e:Event = null):void
		{
			sparkle_canvas.visible = !sparkle_canvas.visible;
		}
		
		protected function goFullScreen(e:Event = null):void
		{
			if (stage.displayState == StageDisplayState.FULL_SCREEN) {
				stage.displayState = StageDisplayState.NORMAL;
			} else {
				try {
					stage.fullScreenSourceRect = new Rectangle(0, 0, 780, 460);
					stage.displayState = 'fullScreen';
				} catch(e:Error){};
			}
		}
		
		protected function onStageKeyDown(e:KeyboardEvent):void
		{
			if(e.keyCode == Keyboard.SPACE)
			{
				toggleSparkles();
			}
		}

		protected function initStage():void
		{
			stage.scaleMode = StageScaleMode.NO_SCALE;
			//stage.align = StageAlign.TOP_LEFT;
			
			stage.addEventListener( KeyboardEvent.KEY_DOWN, onStageKeyDown );
			
			stage.doubleClickEnabled = true;
			stage.addEventListener(MouseEvent.DOUBLE_CLICK, goFullScreen);

			var myContextMenu:ContextMenu = new ContextMenu();
			myContextMenu.hideBuiltInItems();

			var copyr:ContextMenuItem = new ContextMenuItem("© inspirit.ru", true, false);			
			
			var item:ContextMenuItem;
			item = new ContextMenuItem("Load mp3...");
			item.addEventListener(ContextMenuEvent.MENU_ITEM_SELECT, loadMP3, false, 0, true);
			myContextMenu.customItems.push(item);
			
			item = new ContextMenuItem("Switch to auto mode", false, false);
			item.addEventListener(ContextMenuEvent.MENU_ITEM_SELECT, toAutoMode, false, 0, true);
			myContextMenu.customItems.push(item);
			
			item = new ContextMenuItem("Toggle Sparkles");
			item.addEventListener(ContextMenuEvent.MENU_ITEM_SELECT, toggleSparkles, false, 0, true);
			myContextMenu.customItems.push(item);
			
			item = new ContextMenuItem("Fullscreen", true);
			item.addEventListener(ContextMenuEvent.MENU_ITEM_SELECT, goFullScreen, false, 0, true);
			myContextMenu.customItems.push(item);
			
			myContextMenu.customItems.push(copyr);

			contextMenu = myContextMenu;
		}
	}
}
