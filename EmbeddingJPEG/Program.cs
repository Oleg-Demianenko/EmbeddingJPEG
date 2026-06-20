using System;

namespace EmbeddingJPEG
{
    class Program
    {
        static int Main()
        {
            StegoImageJPEG image = new StegoImageJPEG("image.jpg");
            image.Parse();

            string message = "Almost secret message";

            EmbedderJPEG embedder = new EmbedderJPEG(image);
            embedder.EmbedMessage(message, "embedded_image.jpg", 1);

            StegoImageJPEG stegoImage = new StegoImageJPEG("embedded_image.jpg");
            stegoImage.Parse();

            embedder = new EmbedderJPEG(stegoImage);
            string extractedMessage = embedder.ReadMessageFromImage(stegoImage, message.Length, 1);

            Console.WriteLine(extractedMessage);

            return 0;
        }
    }
}