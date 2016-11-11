package 
{
	import example.away3dlite.ExampleFlocking;
	import example.ribbon.ExampleRibbon;

	import com.bit101.components.Label;
	import com.bit101.components.Panel;
	import com.bit101.components.Style;

	import flash.display.Sprite;
	import flash.display.StageScaleMode;
	import flash.events.Event;
	import flash.geom.Point;
	import flash.ui.ContextMenu;
	import flash.ui.ContextMenuItem;
	import flash.utils.getTimer;

	/**
	 * Steerig Behaviors
	 * Example usage. Uncomment away3dliteExample or ribbonExample methods
	 *
	 * @author Eugene Zatepyakin
	 * @link http://blog.inspirit.ru
	 */

	[SWF(width='800',height='640',frameRate='33',backgroundColor='0xFFFFFF')]

	public class Main extends Sprite
	{
		protected static const ORIGIN:Point = new Point();

		protected var view:Sprite;
		
		protected var away3dliteFlockEx:ExampleFlocking;
		protected var ribbonEx:ExampleRibbon;

		protected var p:Panel;
		protected var _timer:uint;
		protected var _fps:uint;
		protected var _ms:uint;
		protected var _ms_prev:uint;

		public function Main()
		{
			if(stage) init();
			else addEventListener(Event.ADDED_TO_STAGE, init);
		}

		private function init(e:Event = null):void
		{
			removeEventListener(Event.ADDED_TO_STAGE, init);
			initStage();
			iniPanel();

			view = new Sprite();
			view.y = 40;
			addChild(view);
			
			//away3dliteExample();
			ribbonExample();

			addEventListener(Event.ENTER_FRAME, render);
		}

		protected function render( e:Event = null ) : void
		{
			//away3dliteFlockEx.render();
			ribbonEx.render();
			
			//Label(p.getChildByName('stat_txt')).text = 'STEER: ' + away3dliteFlockEx.steerTime + 'ms\nRENDER: ' + away3dliteFlockEx.renderTime + 'ms';
			Label(p.getChildByName('stat_txt')).text = 'STEER & RENDER: ' + ribbonEx.steerTime + 'ms';
			
			countFrameTime();
		}
		
		protected function ribbonExample():void
		{
			ribbonEx = new ExampleRibbon();
			view.addChild(ribbonEx);
			//ribbonEx.initFollowers(20);
			ribbonEx.initFlockingBoids(30);
		}

		protected function away3dliteExample():void
		{
			away3dliteFlockEx = new ExampleFlocking();
			view.addChild(away3dliteFlockEx);
			away3dliteFlockEx.initFollowers();
			//away3dliteFlockEx.initFlockingBoids();
		}

		protected function iniPanel():void
		{
			p = new Panel(this);
			p.width = 800;
			p.height = 40;

			Style.PANEL = 0x181f20;
			Style.BUTTON_FACE = 0x181f20;
			Style.LABEL_TEXT = 0xFFFFFF;

			var lb:Label = new Label(p, 10, 5);
			lb.name = 'fps_txt';
			lb = new Label(p, 10, 20);
			lb.name = 'ms_txt';
			lb = new Label(p, 100, 5);
			lb.name = 'stat_txt';
			
		}
		
		protected function countFrameTime(e:Event = null):void
		{
			_timer = getTimer();
			var lab:Label = Label(p.getChildByName('fps_txt'));
			if( _timer - 1000 >= _ms_prev )
			{
				_ms_prev = _timer;

				lab.text = 'FPS: ' + _fps + ' / ' + stage.frameRate;

				_fps = 0;
			}
			
			Label(p.getChildByName('ms_txt')).text = "MS: " + (_timer - _ms);

			_fps ++;
			_ms = _timer;
		}

		protected function initStage():void
		{
			stage.scaleMode = StageScaleMode.NO_SCALE;
			//stage.align = StageAlign.TOP_LEFT;

			var myContextMenu:ContextMenu = new ContextMenu();
			myContextMenu.hideBuiltInItems();


			var copyr:ContextMenuItem = new ContextMenuItem("© inspirit.ru", true, false);
			myContextMenu.customItems.push(copyr);

			contextMenu = myContextMenu;
		}
	}
}
