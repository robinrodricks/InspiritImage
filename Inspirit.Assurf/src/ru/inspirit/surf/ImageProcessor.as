package ru.inspirit.surf 
{
	import flash.geom.Rectangle;
	import flash.display.BitmapData;
	import flash.filters.ColorMatrixFilter;
	import flash.geom.Point;

	/**
	 * Image processor is used to optimize image quality before detection
	 * please note that the result output image should be Grayscale
	 * 
	 * @author Eugene Zatepyakin
	 */
	public class ImageProcessor 
	{
		public static const GRAYSCALE_MATRIX:ColorMatrixFilter = new ColorMatrixFilter([
			.2989, .587, .114, 0, 0,
            .2989, .587, .114, 0, 0,
            .2989, .587, .114, 0, 0,
            0, 0, 0, 0, 0
		]);
		
		public static const ORIGIN:Point = new Point();
		
		public var imageRect:Rectangle;
		
		/**
		 * Pre-process image before detection
		 * 
		 * @param input		source image
		 * @param output	result image
		 */
		public function preProcess(input:BitmapData, output:BitmapData):void
		{
			output.applyFilter(input, imageRect, ORIGIN, GRAYSCALE_MATRIX);
		}
	}
}
