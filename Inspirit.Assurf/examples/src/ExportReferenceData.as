package
{
	import com.bit101.components.RadioButton;
	import ru.inspirit.surf.ASSURF;
	import ru.inspirit.surf.Utils;

	import com.bit101.components.CheckBox;
	import com.bit101.components.HUISlider;
	import com.bit101.components.Label;
	import com.bit101.components.PushButton;
	import com.bit101.components.Text;

	import flash.display.BitmapData;
	import flash.display.LoaderInfo;
	import flash.display.Sprite;
	import flash.display.StageScaleMode;
	import flash.events.Event;
	import flash.utils.ByteArray;
	import flash.utils.Endian;

	/**
	 * @author Eugene Zatepyakin
	 */
	 
	[SWF(width='410',height='200',frameRate='30',backgroundColor='0xFFFFFF')]
	
	public final class ExportReferenceData extends Sprite
	{
		protected var myview:Sprite;
		protected var info:Text;
		
		public const surf:ASSURF = new ASSURF();

		public var detectPrecision:int = ASSURF.DETECT_PRECISION_MEDIUM;
		public var scaleLevels:int = 4;
		public var maxPointsPerLevel:int = 1500;
		public var supressNeighb:Boolean = false;
		public var bmp:BitmapData;
		
		public function ExportReferenceData()
		{
			if(stage) init();
			else addEventListener(Event.ADDED_TO_STAGE, init);
		}
		
		protected function init():void
		{
			removeEventListener(Event.ADDED_TO_STAGE, init);
			
			stage.scaleMode = StageScaleMode.NO_SCALE;

			myview = new Sprite();
			
			var pb:PushButton;
			var lb:Label = new Label( myview, 10, 10, 'SCALE LEVELS' );
			var sc:HUISlider = new HUISlider( myview , 3, 25, '', onScaleLevelChange);
			sc.width = 140;
			sc.setSliderParams(2, 5, scaleLevels);
			sc.labelPrecision = 0;
			
			lb = new Label( myview, 117, 10, 'POINTS PER LEVEL' );
			sc = new HUISlider( myview , 110, 25, '', onPointsPerLevelChange);
			sc.width = 170;
			sc.setSliderParams(500, 2000, maxPointsPerLevel);
			sc.labelPrecision = 0;

			lb = new Label( myview, 270, 10, 'DETECT PRECISION' );
			var chk:RadioButton = new RadioButton( myview , 270, 29, 'HIGH', false, onPrecisionChange);
			chk.groupName = 'precision';
			chk = new RadioButton( myview , 310, 29, 'MEDIUM', true, onPrecisionChange);
			chk.groupName = 'precision';
			chk = new RadioButton( myview , 365, 29, 'LOW', false, onPrecisionChange);
			chk.groupName = 'precision';
			
			new CheckBox( myview , 270, 49, 'SUPRESS NEIGHBORS', onSupressNeighbChange);

			info = new Text( myview , 12, 45, ''); 
			info.height = 50;
			
			pb = new PushButton( myview, 12, 105, 'SELECT IMAGE', onImageSelect);
			pb.width = 100;
			
			pb = new PushButton( myview, 122, 105, 'PROCESS', processImage);
			pb.width = 100;
			
			pb = new PushButton( myview, 232, 105, 'EXPORT DATA', exportData);
			pb.width = 100;
			
			addChild(myview);
		}
		
		protected function onImageSelect(e:Event):void
		{
			Utils.openImageFile(onImageLoaded);
		}
		
		protected function onScaleLevelChange(e:Event):void
		{
			scaleLevels = HUISlider(e.currentTarget).value;
		}
		
		protected function onPointsPerLevelChange(e:Event):void
		{
			maxPointsPerLevel = HUISlider(e.currentTarget).value;
		}
		
		protected function onSupressNeighbChange(e:Event):void
		{
			supressNeighb = CheckBox(e.currentTarget).selected;
		}
		
		protected function onPrecisionChange(e:Event):void
		{
			var chk:RadioButton = RadioButton( e.currentTarget );
			switch(chk.label)
			{
				case 'HIGH':
					detectPrecision = 0;
					break;
				case 'MEDIUM':
					detectPrecision = 1;
					break;
				case 'LOW':
					detectPrecision = 2;
					break;
			}
		}
		
		protected function onImageLoaded(e:Event):void
		{
			var ld:LoaderInfo = LoaderInfo(e.currentTarget);
			ld.removeEventListener(Event.COMPLETE, onImageLoaded);
			
			if(bmp) bmp.dispose();
			bmp = new BitmapData(ld.width, ld.height, false, 0x00);
			bmp.draw(ld.content);
			
			info.text = 'IMAGE LOADED [' + bmp.width + 'x' + bmp.height + ']';
		}
		
		protected function processImage(e:Event):void
		{
			surf.init(detectPrecision, 300, 10000, 1);
			surf.addRefObject( bmp, scaleLevels, maxPointsPerLevel, supressNeighb );
			surf.buildRefIndex();
			
			info.text = 'IMAGE LOADED [' + bmp.width + 'x' + bmp.height + ']';
			info.text += '\n' + surf.referencePointsCount + ' POINTS PROCESSED';
		}
		
		protected function exportData(e:Event):void
		{
			var ba:ByteArray = new ByteArray();
			ba.endian = Endian.LITTLE_ENDIAN;
			surf.exportReferenceData(ba);
			ba.position = 0;
			ba.compress();
			
			info.text = 'IMAGE LOADED [' + bmp.width + 'x' + bmp.height + ']';
			info.text += '\n' + surf.referencePointsCount + ' POINTS PROCESSED';
			info.text += '\n' + 'SAVED ' + int(ba.length/1024+0.5) + 'KB';
			
			Utils.savePointsData(ba);
		}
	}
}
