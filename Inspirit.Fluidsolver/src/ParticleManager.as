package  
{
	import ru.inspirit.utils.ColorUtils;
	import ru.inspirit.utils.NumberUtils;
	import ru.inspirit.utils.Random;

	import flash.display.BitmapData;

	/**
	 * Particle manager implementing simple LinkedList
	 * @author Eugene Zatepyakin
	 */
	public class ParticleManager 
	{
		public static const MAX_PARTICLES:uint = 5000;
		public static const VMAX:Number = 0.013;
		public static const VMAX2:Number = VMAX * VMAX;
		
		
		public var head:Particle;
		public var tail:Particle;
		
		public function ParticleManager() 
		{
			reset();
		}
		
		public function update(bmp:BitmapData, lines:Boolean):void
		{
			var p:Particle = head;			
			
			var a:int;
			var c:uint;
			
			var vxNorm:Number;
			var vyNorm:Number;
			var satInc:Number;
			var v2:Number;
			var m:Number;
			var rgb:Object;
			
			while (p)
			{
				if (p.alpha > 0) 
				{
					p.update();
					
					a = int(p.alpha * 0xFF + .5);
					
					if(Main.drawFluid){
						c = a << 24 | a << 16 | a << 8 | a;
					} else {
						vxNorm = p.vx * Main.isw;
		                vyNorm = p.vy * Main.ish;
		                v2 = vxNorm * vxNorm + vyNorm * vyNorm;
						
						if(v2 > VMAX2) v2 = VMAX2;
						
						m = p.mass;
						satInc = m > 0.5 ? m * m * m : 0;
						satInc *= satInc * satInc * satInc;
						
						rgb = ColorUtils.HSB2GRB(0, NumberUtils.map(v2, 0, VMAX2, 0, 1) + satInc, NumberUtils.interpolate(m, 0.5, 1) * p.alpha);
						c = a << 24 | (rgb.r * 0xFF) << 16 | (rgb.g * 0xFF) << 8 | rgb.b * 0xFF;
					}
					
					line(bmp, int(p.x - p.vx - Main.mx + .5), int(p.y - p.vy - Main.my + .5), int(p.x + .5), int(p.y + .5), c);
					
				}
				p = p.next;
			}
		}
		
		public function reset():void
		{
			var p:Particle = new Particle();
			
			head = tail = p;
			
			var k:int = MAX_PARTICLES;
			
			for (var i:int = 1; i < k; i++)
			{
				p = new Particle();
				tail.next = p;
				tail = p;
			}
		}
		
		public function addParticles(x:Number, y:Number, n:int):void
		{
			while ( --n > -1) {
				addParticle(x + Random.float(-15, 15), y + Random.float(-15, 15));
			}
		}
		
		public function addParticle(x:Number, y:Number):void
		{
			var p:Particle = head;
			
			head = head.next;
			tail.next = p;
			p.next = null;
			tail = p;
			
			tail.init(x, y);
		}
		
		public static function line(bmp:BitmapData, x0:int, y0:int, x1:int, y1:int, c:uint):void
		{	
			var x:int = x0;
			var y:int = y0;		
			var dx:int = x1 - x0;
			var dy:int = y1 - y0;
			var xinc:int = ( dx > 0 ) ? 1 : -1;
			var yinc:int = ( dy > 0 ) ? 1 : -1;
			var cumul:int;
			var i:int;			
			
			dx = (dx ^ (dx >> 31)) - (dx >> 31);//abs
			dy = (dy ^ (dy >> 31)) - (dy >> 31);
			
			bmp.setPixel32(x, y, c);
			
			if ( dx > dy ) {
				cumul = dx >> 1 ;
		  		for ( i = 1; i <= dx; ++i ) {
					x += xinc;
					cumul += dy;
					if (cumul >= dx){
			  			cumul -= dx;
			  			y += yinc;
					}
					bmp.setPixel32(x, y, c);
				}
			} else {
		  		cumul = dy >> 1;
		  		for ( i = 1; i <= dy; ++i ) {
					y += yinc;
					cumul += dx;
					if ( cumul >= dy ) {
			  			cumul -= dy;
			  			x += xinc;
					}
					bmp.setPixel32(x, y, c);
				}
			}
		}
		
	}
	
}