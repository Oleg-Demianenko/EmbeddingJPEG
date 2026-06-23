using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace EmbeddingJPEG
{
    public class EmbedderJPEG
    {
        StegoImageJPEG image;
        static int[,] H1b;
        static int[,] H2b;
        static int[] H1;
        static int[] H2;
        Dictionary<int, int> HCols;
        Dictionary<int, int> BCHCols;
        Dictionary<int, int> BCHColsComb2;
        Dictionary<int, int> BCHColsComb3;
        int maxBlockCount;
        public EmbedderJPEG(StegoImageJPEG image)
        {
            this.image = image;
            HCols = new Dictionary<int, int>();
            BCHCols = new Dictionary<int, int>();
            BCHColsComb2 = new Dictionary<int, int>();
            BCHColsComb3 = new Dictionary<int, int>();
            MakeMatrices();
        }
        /// <summary>
        /// Представление проверочных матриц кодов в удобном для работы виде
        /// </summary>
        private void MakeMatrices()
        {
            H1b = new int[,]
            {
            { 1, 0, 0, 0, 0, 0, 0, 0 },
            { 1, 1, 0, 0, 0, 0, 0, 0 },
            { 0, 1, 1, 0, 0, 0, 0, 0 },
            { 1, 0, 1, 1, 0, 0, 0, 0 },
            { 0, 1, 0, 1, 1, 0, 0, 0 },
            { 0, 0, 1, 0, 1, 1, 0, 0 },
            { 0, 0, 0, 1, 0, 1, 1, 0 },
            { 1, 0, 0, 0, 1, 0, 1, 1 },
            { 0, 1, 0, 0, 0, 1, 0, 1 },
            { 0, 0, 1, 0, 0, 0, 1, 0 },
            { 0, 0, 0, 1, 0, 0, 0, 1 },
            { 0, 0, 0, 0, 1, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 1, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 1, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 1 },
            };
            H2b = new int[,]
            {
            { 1, 0, 0 },
            { 0, 1, 0 },
            { 1, 1, 0 },
            { 0, 0, 1 },
            { 1, 0, 1 },
            { 0, 1, 1 },
            { 1, 1, 1 }
            };
            H1 = new int[8];
            H2 = new int[3];
            // Записываем столбцы транспонированных матриц как числа
            int sum;
            for (int j = 0; j < H1b.GetLength(1); j++)
            {
                sum = 0;
                for (int i = 0; i < H1b.GetLength(0); i++)
                    sum = (sum << 1) | H1b[i, j];
                H1[j] = sum;
            }
            for (int j = 0; j < H2b.GetLength(1); j++)
            {
                sum = 0;
                for (int i = 0; i < H2b.GetLength(0); i++)
                    sum = (sum << 1) | H2b[i, j];
                H2[j] = sum;
            }
            // Записываем столбцы не транспонированных матриц как числа
            for (int i = 0; i < H2b.GetLength(0); i++)
            {
                sum = 0;
                for (int j = 0; j < H2b.GetLength(1); j++)
                    sum = (sum << 1) | H2b[i, j];
                HCols.Add(sum, i);
            }
            for (int i = 0; i < H1b.GetLength(0); i++)
            {
                sum = 0;
                for (int j = 0; j < H1b.GetLength(1); j++)
                    sum = (sum << 1) | H1b[i, j];
                BCHCols.Add(sum, i);
            }
            // Собираем комбинации длины 2 и 3 для БЧХ
            foreach (int v1 in BCHCols.Keys)
                foreach (int v2 in BCHCols.Keys)
                {
                    sum = v1 ^ v2;
                    if (sum != 0 && !BCHCols.ContainsKey(sum) && !BCHColsComb2.ContainsKey(sum))
                        BCHColsComb2.Add(sum, (BCHCols[v1] << 4) | BCHCols[v2]);
                }
            foreach (int v1 in BCHCols.Keys)
                foreach (int v2 in BCHColsComb2.Keys)
                {
                    sum = v1 ^ v2;
                    if (sum != 0 && !BCHCols.ContainsKey(sum) &&
                        !BCHColsComb2.ContainsKey(sum) && !BCHColsComb3.ContainsKey(sum))
                        BCHColsComb3.Add(sum, (BCHCols[v1] << 8) | BCHColsComb2[v2]);
                }
        }
        /// <summary>
        /// Выбор коэффициентов, находящихся в диапазоне рабочей области
        /// </summary>
        /// <param name="DCTCoeffs">Массив коэффициентов ДКП</param>
        /// <param name="blockCount">Число блоков, которые требуется набрать</param>
        /// <returns>Массив коэффициентов, разледённых на две зоны</returns>
        private int[,] SelectDCTCoeffs(short[,,] DCTCoeffs, int blockCount)
        {
            int[,] coeffs = new int[4, blockCount * 22];
            int blockIndex = 0;
            int coeffIndex = 0;
            int arrayCoeffIndex = 0;
            int coeff;
            for (int m = 0; m < DCTCoeffs.GetLength(0); m++)
                for (int b = 0; b < DCTCoeffs.GetLength(1); b++)
                    for (int c = 1; c < DCTCoeffs.GetLength(2); c++)
                    {
                        coeff = Math.Abs(DCTCoeffs[m, b, c]);
                        if (coeff == 1) // {-1; 1}
                        {
                            arrayCoeffIndex = blockIndex * 22 + coeffIndex;
                            coeffs[0, arrayCoeffIndex] = DCTCoeffs[m, b, c];
                            coeffs[1, arrayCoeffIndex] = m;
                            coeffs[2, arrayCoeffIndex] = b;
                            coeffs[3, arrayCoeffIndex] = c;
                            coeffIndex++;
                        }
                        if (coeffIndex == 15)
                        {
                            blockIndex++;
                            coeffIndex = 0;
                        }
                        if (blockIndex == blockCount)
                            goto EndLightBlocks;
                    }
                EndLightBlocks:;
            coeffIndex = 15;
            blockIndex = 0;
            for (int m = 0; m < DCTCoeffs.GetLength(0); m++)
                for (int b = 0; b < DCTCoeffs.GetLength(1); b++)
                    for (int c = 1; c < DCTCoeffs.GetLength(2); c++)
                    {
                        coeff = Math.Abs(DCTCoeffs[m, b, c]);
                        if (coeff > 1 && coeff < 8) // [-7; -2] + [2; 7]
                        {
                            arrayCoeffIndex = blockIndex * 22 + coeffIndex;
                            coeffs[0, arrayCoeffIndex] = DCTCoeffs[m, b, c];
                            coeffs[1, arrayCoeffIndex] = m;
                            coeffs[2, arrayCoeffIndex] = b;
                            coeffs[3, arrayCoeffIndex] = c;
                            coeffIndex++;
                        }
                        if (coeffIndex == 22)
                        {
                            blockIndex++;
                            coeffIndex = 15;
                        }
                        if (blockIndex == blockCount)
                            goto EndHeavyBlocks;
                    }
                EndHeavyBlocks:;
            return coeffs;
        }
        /// <summary>
        /// Метод, собирающий блоки контейнера из младших бит коэффициентов
        /// </summary>
        /// <param name="coeffs">Массив коэффициентов рабочей области</param>
        /// <param name="blockCount">Количество блоков контейнера</param>
        /// <returns>Массив блоков контейнера, разделённый на две зоны</returns>
        private int[,] CollectContainerBlocks(int[,] coeffs, int blockCount)
        {
            int[,] blocks = new int[2, blockCount];
            int sum = 0;
            int coeffIndex = 0;
            for (int b = 0; b < blockCount; b++)
            {
                sum = 0;
                for (int k = 0; k < 15; k++)
                {   // -1 -> 0; 1 -> 1
                    sum = (sum << 1) | (coeffs[0, coeffIndex] == 1 ? 1 : 0);
                    coeffIndex++;
                }
                blocks[0, b] = sum;
                sum = 0;
                for (int k = 0; k < 7; k++)
                {
                    sum = (sum << 1) | (Math.Abs(coeffs[0, coeffIndex]) & 1);
                    coeffIndex++;
                }
                blocks[1, b] = sum;
            }
            return blocks;
        }
        /// <summary>
        /// Метод, собирающий блоки контейнера из младших бит коэффициентов.
        /// Набирает последовательности блоков, распределённые по рабочей области,
        /// в количестве равном числу копий сообщения
        /// </summary>
        /// <param name="coeffs">Массив коэффициентов рабочей области</param>
        /// <param name="blockCount">Количество блоков контейнера</param>
        /// <returns>Массив блоков контейнера, разделённый на две зоны</returns>
        private int[,] CollectContainerBlocks(int[,] coeffs, int blockCount, int copiesCount)
        {
            int[,] blocks = new int[2, blockCount * copiesCount];
            int gapBase = 0;
            int gapRem = 0;
            if (copiesCount > 1)
            {
                gapBase = (maxBlockCount - blockCount * copiesCount) / (copiesCount - 1);
                gapRem = (maxBlockCount - blockCount * copiesCount) % (copiesCount - 1);
            }
            int coeffIndex = 0;
            int sum;
            for (int cop = 0; cop < copiesCount; cop++)
            {
                for (int b = 0; b < blockCount; b++)
                {
                    sum = 0;
                    for (int k = 0; k < 15; k++)
                    {
                        sum = (sum << 1) | (coeffs[0, coeffIndex] == 1 ? 1 : 0);
                        coeffIndex++;
                    }
                    blocks[0, cop * blockCount + b] = sum;
                    sum = 0;
                    for (int k = 0; k < 7; k++)
                    {
                        sum = (sum << 1) | (Math.Abs(coeffs[0, coeffIndex]) & 1);
                        coeffIndex++;
                    }
                    blocks[1, cop * blockCount + b] = sum;
                }
                if (gapRem > 0)
                    coeffIndex += (gapBase + 1) * 22;
                else
                    coeffIndex += gapBase * 22;
                gapRem--;
            }
            return blocks;
        }
        /// <summary>
        /// Вспомогательный метод для переводя строки сообщения в байты
        /// </summary>
        /// <param name="message">Сообщение</param>
        /// <returns>Байты сообщения</returns>
        private byte[] MessageBytes(string message)
        {
            return Encoding.ASCII.GetBytes(message);
        }
        /// <summary>
        /// Добавление байт валидации
        /// </summary>
        /// <param name="messageBytes">Массив байт сообщения (минимум два свободных)</param>
        /// <param name="length">Длина сообщения в байтах</param>
        private void AddCRC(byte[] messageBytes, int length)
        {
            byte[] buffer = messageBytes.ToArray();
            int crcPolynom = 0b_10001000000100001;
            int restBytes = length;
            int m = (messageBytes[0] << 8) | messageBytes[1];
            int window = 0;
            int byteIndex = 2;
            while (restBytes > 0)
            {
                window = 0;
                m = ((m << 8) & 0x00FFFFFF) | buffer[byteIndex];
                for (int i = 0; i < 8; i++)
                {
                    window = m & (0b_11111111111111111 << (7 - i));
                    if ((window >> (23 - i)) == 1)
                        m = m ^ (crcPolynom << (7 - i));
                }
                buffer[byteIndex - 2] = (byte)(m >> 16);
                byteIndex++;
                restBytes--;
            }
            messageBytes[byteIndex - 2] = (byte)((m >> 8) & 0xFF);
            messageBytes[byteIndex - 1] = (byte)(m & 0xFF);
        }
        /// <summary>
        /// Добавление маркеров сихронизации
        /// </summary>
        /// <param name="messageBytes">Массив байт сообщения (минимум 4 свободных)</param>
        /// <param name="index">Индекс копии сообщения</param>
        private void AddSyncMarker(byte[] messageBytes, int index)
        {
            messageBytes[messageBytes.Length - 4] = 0xF9;
            messageBytes[messageBytes.Length - 3] = (byte)(0xA8 + index);
            messageBytes[messageBytes.Length - 2] = 0xE0;
            messageBytes[messageBytes.Length - 1] = 0x8A;
        }
        /// <summary>
        /// Метод, набирающий блоки сообщения из байт сообщения. Добавляет
        /// байты валидации и маркеры синхронизации
        /// </summary>
        /// <param name="messageBytes">Байты сообщения</param>
        /// <param name="blockCount">Количество блоков</param>
        /// <param name="copiesCount">Количество копий сообщения</param>
        /// <returns></returns>
        private int[,] GetMessageBlocks(byte[] messageBytes, int blockCount, int copiesCount)
        {
            byte[] bytes = messageBytes.Concat(new byte[6]).ToArray();
            AddCRC(bytes, bytes.Length - 6);
            int[,] blocks = new int[2, blockCount * copiesCount];
            int getMask = 0;
            int setMask = 0;
            int buffer = 0;
            int length = 0;
            int byteIndex = 0;
            int blockIndex = 0;
            for (int cop = 0; cop < copiesCount; cop++)
            {
                AddSyncMarker(bytes, cop);
                byteIndex = 0;
                for (int i = 0; i < blockCount - 4; i++)
                {
                    if (length == 8)
                    {
                        blocks[0, blockIndex] = buffer;
                        buffer = 0;
                        length = 0;
                    }
                    else
                    {
                        if (byteIndex < bytes.Count())
                        {
                            buffer = (buffer << 8) | bytes[byteIndex];
                            byteIndex++;
                            getMask = 0b_11111111 << (length);
                            setMask = (1 << (length)) - 1;
                            blocks[0, blockIndex] = (buffer & getMask) >> length;
                            buffer = buffer & setMask;
                        }
                        else
                        {
                            blocks[0, blockIndex] = buffer << (8 - length);
                            buffer = 0;
                            length = 0;
                        }
                    }
                    if (length >= 3)
                    {
                        getMask = 0b_111 << (length - 3);
                        setMask = (1 << (length - 3)) - 1;
                        blocks[1, blockIndex] = (buffer & getMask) >> (length - 3);
                        buffer = buffer & setMask;
                        length -= 3;
                    }
                    else
                    {
                        if (byteIndex < bytes.Count())
                        {
                            buffer = (buffer << 8) | bytes[byteIndex];
                            byteIndex++;
                            getMask = 0b_111 << (length + 5);
                            setMask = (1 << (length + 5)) - 1;
                            blocks[1, blockIndex] = (buffer & getMask) >> (length + 5);
                            buffer = buffer & setMask;
                            length += 5;
                        }
                        else
                        {
                            blocks[1, blockIndex] = buffer << (3 - length);
                        }
                    }
                    blockIndex++;
                }
                // Первая зона - два байта валидации, два байта маркера синхронизации
                // Вторая зона - 4 блока по три бита маркера синхронизации
                length = 0;
                blocks[0, blockIndex] = bytes[bytes.Length - 6];
                blocks[1, blockIndex] = (bytes[bytes.Length - 2] >> 5) & 7;
                blockIndex++;
                blocks[0, blockIndex] = bytes[bytes.Length - 5];
                blocks[1, blockIndex] = (bytes[bytes.Length - 2] >> 1) & 7;
                blockIndex++;
                blocks[0, blockIndex] = bytes[bytes.Length - 4];
                blocks[1, blockIndex] = (bytes[bytes.Length - 1] >> 5) & 7;
                blockIndex++;
                blocks[0, blockIndex] = bytes[bytes.Length - 3];
                blocks[1, blockIndex] = (bytes[bytes.Length - 1] >> 1) & 7;
                blockIndex++;
            }
            return blocks;
        }
        /// <summary>
        /// Синдромное встраивание блоков сообщения в блоки контейнера
        /// </summary>
        /// <param name="containerBlocks">Массив блоков контейнера</param>
        /// <param name="messageBlocks">Массив блоков сообщения</param>
        /// <returns>Массив заполенных блоков контейнера</returns>
        private int[,] EmbedSyndromes(int[,] containerBlocks, int[,] messageBlocks)
        {
            int[,] filledBlocks = new int[containerBlocks.GetLength(0), containerBlocks.GetLength(1)];
            int blockCount = messageBlocks.GetLength(1);
            int contSyndrome;
            int errSyndrome;
            int error;
            for (int i = 0; i < blockCount; i++)
            {
                contSyndrome = 0;
                for (int j = 0; j < H1.Count(); j++)
                    contSyndrome = (contSyndrome << 1) | (BitOperations.PopCount((uint)(containerBlocks[0, i] & H1[j])) & 1);
                errSyndrome = contSyndrome ^ messageBlocks[0, i];
                if (errSyndrome == 0) // Синдром уже совпадает с сообщением
                    error = 0;
                else if (BCHCols.ContainsKey(errSyndrome)) // Достаточно одной ошибки
                    error = 1 << (14 - BCHCols[errSyndrome]);
                else if (BCHColsComb2.ContainsKey(errSyndrome)) // Достаточно двух ошибок
                    error = (1 << (14 - (BCHColsComb2[errSyndrome] >> 4))) |
                        (1 << (14 - (BCHColsComb2[errSyndrome] & 0b_00001111)));
                else // Необходимо добавить три ошибки к блоку контейнера
                    error = (1 << (14 - (BCHColsComb3[errSyndrome] >> 8))) |
                        (1 << (14 - ((BCHColsComb3[errSyndrome] >> 4) & 0b_00001111))) |
                        (1 << (14 - (BCHColsComb3[errSyndrome] & 0b_000000001111)));
                filledBlocks[0, i] = containerBlocks[0, i] ^ error;

                contSyndrome = 0;
                for (int j = 0; j < H2.Count(); j++)
                    contSyndrome = (contSyndrome << 1) | (BitOperations.PopCount((uint)(containerBlocks[1, i] & H2[j])) & 1);
                errSyndrome = contSyndrome ^ messageBlocks[1, i];
                if (errSyndrome == 0) // Синдром уже совпадает с сообщением
                    error = 0;
                else // Необходимо добавить одну ошибку к блоку контейнера
                    error = 1 << (6 - HCols[errSyndrome]);
                filledBlocks[1, i] = containerBlocks[1, i] ^ error;
            }
            for (int i = blockCount; i < containerBlocks.GetLength(1); i++)
            {
                filledBlocks[0, i] = containerBlocks[0, i];
                filledBlocks[1, i] = containerBlocks[1, i];
            }
            return filledBlocks;
        }
        /// <summary>
        /// Разбор блоков контейнера на коэффициенты
        /// </summary>
        /// <param name="containerBlocks">Массив блоков контейнера</param>
        /// <param name="oldCoeffs">Массив старых коэффициентов</param>
        /// <param name="blockCount">Количество блоков</param>
        /// <param name="copiesCount">Количество копий сообщения</param>
        /// <returns>Массив коэффициентов, содержащих биты контейнера</returns>
        int[,] CoeffsFromContainerBlocks(int[,] containerBlocks, int[,] oldCoeffs, int blockCount, int copiesCount)
        {
            int[,] newCoeffs = new int[oldCoeffs.GetLength(0), oldCoeffs.GetLength(1)];
            int gapBase = 0;
            int gapRem = 0;
            if (copiesCount > 1)
            {
                gapBase = (maxBlockCount - blockCount * copiesCount) / (copiesCount - 1);
                gapRem = (maxBlockCount - blockCount * copiesCount) % (copiesCount - 1);
            }
            int contBlock;
            int coeffIndex = 0;
            int targetIndex = 0;

            for (int cop = 0; cop < copiesCount; cop++)
            {
                for (int b = 0; b < blockCount; b++)
                {
                    contBlock = containerBlocks[0, cop * blockCount + b];
                    for (int c = 14; c >= 0; c--)
                    {
                        newCoeffs[0, coeffIndex] = ((contBlock >> c) & 1) == 1 ? 1 : -1;
                        newCoeffs[1, coeffIndex] = oldCoeffs[1, coeffIndex];
                        newCoeffs[2, coeffIndex] = oldCoeffs[2, coeffIndex];
                        newCoeffs[3, coeffIndex] = oldCoeffs[3, coeffIndex];
                        coeffIndex++;
                    }
                    contBlock = containerBlocks[1, cop * blockCount + b];
                    for (int c = 6; c >= 0; c--)
                    {
                        if (oldCoeffs[0, coeffIndex] > 0)
                            newCoeffs[0, coeffIndex] = ((contBlock >> c) & 1) == 1 ? oldCoeffs[0, coeffIndex] | 1 : oldCoeffs[0, coeffIndex] & ~1;
                        else
                            newCoeffs[0, coeffIndex] = ((contBlock >> c) & 1) == 1 ? -((-oldCoeffs[0, coeffIndex]) | 1) : -((-oldCoeffs[0, coeffIndex]) & ~1);
                        newCoeffs[1, coeffIndex] = oldCoeffs[1, coeffIndex];
                        newCoeffs[2, coeffIndex] = oldCoeffs[2, coeffIndex];
                        newCoeffs[3, coeffIndex] = oldCoeffs[3, coeffIndex];
                        coeffIndex++;
                    }
                }
                // Промежутки междку копиями сообщения для их распределения по рабочей области
                if (cop < copiesCount - 1)
                    if (gapRem > 0)
                        targetIndex = coeffIndex + (gapBase + 1) * 22;
                    else
                        targetIndex = coeffIndex + gapBase * 22;
                else
                    targetIndex = newCoeffs.GetLength(1);
                gapRem--;
                while (coeffIndex < targetIndex)
                {
                    newCoeffs[0, coeffIndex] = oldCoeffs[0, coeffIndex];
                    newCoeffs[1, coeffIndex] = oldCoeffs[1, coeffIndex];
                    newCoeffs[2, coeffIndex] = oldCoeffs[2, coeffIndex];
                    newCoeffs[3, coeffIndex] = oldCoeffs[3, coeffIndex];
                    coeffIndex++;
                }
            }
            return newCoeffs;
        }
        /// <summary>
        /// Поиск последовательностей блоков, содержащих копии сообщения
        /// </summary>
        /// <param name="stegoblocks">Массив блоков контейнера</param>
        /// <param name="messageLength">Длина сообщения в байтах</param>
        /// <returns>Массив блоков, содержащих копии сообщения</returns>
        private int[,] SearchFilledBlocks(int[,] stegoblocks, int messageLength)
        {
            List<(int, int)> markers0 = new List<(int, int)>();
            List<(int, int)> markers1 = new List<(int, int)>();
            int buffer1 = 0;
            int buffer2 = 0;
            int extracted = 0;
            int markerIndex = 0;
            int marker0 = 0xF9A8;
            int marker1 = 0xE08A;
            int length = (int)Math.Ceiling(messageLength * 8.0 / 11.0) + 4;

            int counter1 = 0;
            int counter2 = 0;

            for (int b = 0; b < stegoblocks.GetLength(1) - 3; b++)
            {
                for (int i = 0; i < 15; i++)
                {
                    buffer1 = stegoblocks[0, b];
                    buffer1 = buffer1 & ((1 << (15 - i)) - 1);
                    buffer1 = buffer1 << (15 + i);
                    buffer1 = buffer1 | (stegoblocks[0, b + 1] << i);
                    buffer1 = buffer1 | (stegoblocks[0, b + 2] >> (15 - i));
                    buffer2 = buffer1 & ((1 << 15) - 1);
                    buffer1 = buffer1 >> 15;

                    extracted = 0;
                    for (int j = 0; j < H1.Count(); j++)
                        extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer1) & H1[j])) & 1);
                    for (int j = 0; j < H1.Count(); j++)
                        extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer2) & H1[j])) & 1);

                    if (BitOperations.PopCount((uint)(extracted ^ marker0)) < 1)
                    {
                        if (b >= length - 2)
                        {
                            markers0.Add((b + 2 - length, i));
                            counter1++;
                            marker0++;
                        }

                    }
                }
            }
            for (int b = 0; b < stegoblocks.GetLength(1) - 5; b++)
            {
                for (int i = 0; i < 7; i++)
                {
                    buffer1 = stegoblocks[1, b];
                    buffer1 = buffer1 & ((1 << (7 - i)) - 1);
                    buffer1 = buffer1 << (21 + i);
                    buffer1 = buffer1 | (stegoblocks[1, b + 1] << (14 + i));
                    buffer1 = buffer1 | (stegoblocks[1, b + 2] << (7 + i));
                    buffer1 = buffer1 | (stegoblocks[1, b + 3] << i);
                    buffer1 = buffer1 | (stegoblocks[1, b + 4] >> (7 - i));
                    buffer2 = buffer1 & ((1 << 21) - 1);
                    buffer1 = buffer1 >> 21;
                    extracted = 0;
                    for (int j = 0; j < H2.Count(); j++)
                        extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer1) & H2[j])) & 1);
                    extracted = (extracted << 1);
                    buffer1 = buffer2 >> 14;
                    buffer2 = buffer2 & ((1 << 14) - 1);

                    for (int j = 0; j < H2.Count(); j++)
                        extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer1) & H2[j])) & 1);
                    extracted = (extracted << 1);
                    buffer1 = buffer2 >> 7;
                    buffer2 = buffer2 & ((1 << 7) - 1);

                    for (int j = 0; j < H2.Count(); j++)
                        extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer1) & H2[j])) & 1);
                    extracted = (extracted << 1);
                    for (int j = 0; j < H2.Count(); j++)
                        extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer2) & H2[j])) & 1);
                    extracted = (extracted << 1);
                    if (BitOperations.PopCount((uint)(extracted ^ marker1)) < 1)
                    {
                        if (b >= length - 4)
                        {
                            markers1.Add((b + 4 - length, i));
                            counter2++;
                        }
                    }
                }
            }

            buffer1 = stegoblocks[0, stegoblocks.GetLength(1) - 2];
            buffer2 = stegoblocks[0, stegoblocks.GetLength(1) - 1];
            extracted = 0;
            for (int j = 0; j < H1.Count(); j++)
                extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer1) & H1[j])) & 1);
            for (int j = 0; j < H1.Count(); j++)
                extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer2) & H1[j])) & 1);
            if (BitOperations.PopCount((uint)(extracted ^ marker0)) < 1)
            {
                markers0.Add((stegoblocks.GetLength(1) - length, 0));
                counter1++;
                marker0++;
            }

            buffer1 = stegoblocks[1, stegoblocks.GetLength(1) - 4];
            extracted = 0;
            for (int j = 0; j < H2.Count(); j++)
                extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer1) & H2[j])) & 1);
            extracted = (extracted << 1);
            buffer1 = stegoblocks[1, stegoblocks.GetLength(1) - 3];

            for (int j = 0; j < H2.Count(); j++)
                extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer1) & H2[j])) & 1);
            extracted = (extracted << 1);
            buffer1 = stegoblocks[1, stegoblocks.GetLength(1) - 2];

            for (int j = 0; j < H2.Count(); j++)
                extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer1) & H2[j])) & 1);
            extracted = (extracted << 1);
            buffer1 = stegoblocks[1, stegoblocks.GetLength(1) - 1];

            for (int j = 0; j < H2.Count(); j++)
                extracted = (extracted << 1) | (BitOperations.PopCount((uint)((buffer1) & H2[j])) & 1);
            extracted = (extracted << 1);
            if (BitOperations.PopCount((uint)(extracted ^ marker1)) < 1)
            {
                markers1.Add((stegoblocks.GetLength(1) - length, 0));
                counter2++;
            }

            int copiesCount = Math.Min(markers0.Count(), markers1.Count());
            int[,] markerPairs = new int[4, copiesCount];

            if (markers0.Count() <= markers1.Count())
            {
                for (int i = 0; i < markers0.Count(); i++)
                {
                    (int, int) marker0Pos = markers0[i];
                    markerPairs[0, i] = marker0Pos.Item1;
                    markerPairs[1, i] = marker0Pos.Item2;
                    markerPairs[2, i] = int.MaxValue;
                    for (int j = 0; j < markers1.Count(); j++)
                    {
                        (int, int) marker1Pos = markers1[j];
                        if (Math.Abs(markerPairs[0, i] - marker1Pos.Item1) <
                            Math.Abs(markerPairs[0, i] - markerPairs[2, i]))
                        {
                            markerPairs[2, i] = marker1Pos.Item1;
                            markerPairs[3, i] = marker1Pos.Item2;
                        }
                    }
                }
            }
            else
            {
                for (int i = 0; i < markers1.Count(); i++)
                {
                    (int, int) marker1Pos = markers1[i];
                    markerPairs[2, i] = marker1Pos.Item1;
                    markerPairs[3, i] = marker1Pos.Item2;
                    markerPairs[0, i] = int.MaxValue;

                    for (int j = 0; j < markers0.Count(); j++)
                    {
                        (int, int) marker0Pos = markers0[j];
                        if (Math.Abs(markerPairs[2, i] - marker0Pos.Item1) <
                            Math.Abs(markerPairs[2, i] - markerPairs[0, i]))
                        {
                            markerPairs[0, i] = marker0Pos.Item1;
                            markerPairs[0, i] = marker0Pos.Item2;
                        }
                    }
                }
            }

            int[,] blocks = new int[2, copiesCount * length];
            int blockIndex = 0;
            int currentBlock = 0;
            int startBlock = 0;
            int startCoeff = 0;

            for (int cop = 0; cop < copiesCount; cop++)
            {
                startBlock = markerPairs[0, cop];
                startCoeff = markerPairs[1, cop];
                for (int b = startBlock; b < startBlock + length; b++)
                {
                    currentBlock = (stegoblocks[0, b] & ((1 << (15 - startCoeff)) - 1)) << startCoeff;
                    if (b < stegoblocks.GetLength(1) - 1)
                        currentBlock = currentBlock | (stegoblocks[0, b + 1] >> (15 - startCoeff));
                    blocks[0, blockIndex] = currentBlock;
                    blockIndex++;
                }

                startBlock = markerPairs[2, cop];
                startCoeff = markerPairs[3, cop];
                blockIndex -= length;
                for (int b = startBlock; b < startBlock + length; b++)
                {
                    currentBlock = (stegoblocks[1, b] & ((1 << (7 - startCoeff)) - 1)) << startCoeff;
                    if (b < stegoblocks.GetLength(1) - 1)
                        currentBlock = currentBlock | (stegoblocks[1, b + 1] >> (7 - startCoeff));
                    blocks[1, blockIndex] = currentBlock;
                    blockIndex++;
                }
            }
            return blocks;
        }
        /// <summary>
        /// Проверка валидации копии сообщения
        /// </summary>
        /// <param name="bytes">Массив байт сообщения (с байтами валидации)</param>
        /// <returns>1 - успешно, 0 - неуспешно</returns>
        private int CheckValidation(byte[] bytes)
        {
            int blockIndex = bytes.Length - 2;
            int currentRemainder = (bytes[blockIndex] << 8) | bytes[blockIndex + 1];
            AddCRC(bytes, bytes.Length - 2);
            int expectedRemainder = (bytes[blockIndex] << 8) | bytes[blockIndex + 1];
            if (expectedRemainder == 0)
                return 1;
            else return 0;
        }
        /// <summary>
        /// Извлечение копий сообщения из блоков контейнера
        /// </summary>
        /// <param name="stegoBlocks">Массив блоков контейнера со встроенным сообщением</param>
        /// <param name="messageLength">Длина сообщения в байтах</param>
        /// <param name="copiesCount">Количество копий сообщения</param>
        /// <returns>Массив байт подряд идущих копий сообщения,
        /// Массив результатов проверки валидации каждой копии</returns>
        private (byte[], int[]) ExtractMessageCopies(int[,] stegoBlocks, int messageLength, int copiesCount)
        {
            int copyLength = messageLength + 2;
            byte[] copies = new byte[copyLength * copiesCount];
            int blockCount = (int)Math.Ceiling(messageLength * 8.0 / 11.0) + 4;
            int buffer = 0;
            int byteIndex = 0;
            int blockIndex = 0;
            int length = 0;
            int crcRemainder = 0;
            int[] validationResults = new int[copiesCount];
            for (int cop = 0; cop < copiesCount; cop++)
            {
                blockIndex = blockCount * cop;
                length = 0;
                while (byteIndex < copyLength * cop + messageLength)
                {

                    if (length < 8)
                    {
                        for (int i = 0; i < H1.Count(); i++)
                            buffer = (buffer << 1) | (BitOperations.PopCount((uint)(stegoBlocks[0, blockIndex] & H1[i])) & 1);
                        for (int i = 0; i < H2.Count(); i++)
                            buffer = (buffer << 1) | (BitOperations.PopCount((uint)(stegoBlocks[1, blockIndex] & H2[i])) & 1);
                        length += 11;
                        blockIndex++;
                    }
                    copies[byteIndex] = (byte)(buffer >> (length - 8));
                    buffer = buffer & ((1 << (length - 8)) - 1);
                    length -= 8;
                    byteIndex++;
                }
                buffer = 0;
                for (int i = 0; i < H1.Count(); i++)
                    buffer = (buffer << 1) | (BitOperations.PopCount((uint)(stegoBlocks[0, blockIndex] & H1[i])) & 1);
                copies[byteIndex] = (byte)buffer;
                buffer = 0;
                byteIndex++;
                blockIndex++;
                for (int i = 0; i < H1.Count(); i++)
                    buffer = (buffer << 1) | (BitOperations.PopCount((uint)(stegoBlocks[0, blockIndex] & H1[i])) & 1);
                copies[byteIndex] = (byte)buffer;
                buffer = 0;
                length = 0;
                byteIndex++;
                blockIndex++;
                validationResults[cop] = CheckValidation(copies[(copyLength * cop)..(copyLength * (cop + 1))].ToArray());
            }
            return (copies, validationResults);
        }
        /// <summary>
        /// Декодирование сообщения из нескольких копий методом голосвания
        /// </summary>
        /// <param name="copies">Массив байт подряд идущих копий сообщения</param>
        /// <param name="validation">Массив результатов валидации копий</param>
        /// <param name="messageLength">Длина сообщения в байтах</param>
        /// <param name="copiesCount">Количество копий</param>
        /// <returns>Байты сообщения</returns>
        private byte[] DecodeMessage(byte[] copies, int[] validation, int messageLength, int copiesCount)
        {
            byte[] messageBytes = new byte[messageLength];
            int buffer = 0;
            int sum = 0;
            for (int c = 0; c < copiesCount; c++)
                if (validation[c] == 1) // Если нашлась валидная копия - берём её
                {
                    for (int i = 0; i < messageLength; i++)
                    {
                        messageBytes[i] = copies[c * (messageLength + 2) + i];
                    }
                    return messageBytes;
                }
            // Если не нашлось валидной копии, проводим голосование
            for (int i = 0; i < messageLength; i++)
            {
                buffer = 0;
                for (int j = 0; j < 8; j++)
                {
                    sum = 0;
                    for (int c = 0; c < copiesCount; c++)
                        sum += (copies[(messageLength + 2) * c + i] >> (7 - j)) & 1;
                    buffer = buffer << 1;
                    if (sum > copiesCount / 2)
                        buffer++;
                }
                messageBytes[i] = (byte)buffer;
            }
            return messageBytes;
        }
        /// <summary>
        /// Запись новых коэффициентов в общий массив
        /// </summary>
        /// <param name="oldCoeffs">Массив старых коэффициентов рабочей области</param>
        /// <param name="newCoeffs">Массив новых коэффициентов рабочей области</param>
        /// <returns>Новый массив коэффициентов ДКП</returns>
        private short[,,] ReplaceDCTCoeffs(short[,,] oldCoeffs, int[,] newCoeffs)
        {
            short[,,] newDCTCoeff = new short[oldCoeffs.GetLength(0), oldCoeffs.GetLength(1), oldCoeffs.GetLength(2)];
            Array.Copy(oldCoeffs, newDCTCoeff, oldCoeffs.Length);
            for (int i = 0; i < newCoeffs.GetLength(1); i++)
                newDCTCoeff[newCoeffs[1, i], newCoeffs[2, i], newCoeffs[3, i]] = (short)newCoeffs[0, i];
            return newDCTCoeff;
        }
        /// <summary>
        /// Извлечение сообщения из изображения
        /// </summary>
        /// <param name="stegoImage">Изображениесо встроенным сообщением</param>
        /// <param name="length">Длина сообщения в байтах</param>
        /// <param name="n">Количество копий сообщения</param>
        /// <returns>Строка сообщения</returns>
        public string ReadMessageFromImage(StegoImageJPEG stegoImage, int length, int n)
        {
            short[,,] DCTCoeffs = stegoImage.DCTCoeffs;
            int blockCount = (int)Math.Ceiling(length * 8.0 / 11.0) + 4;

            int counter1 = 0;
            int counter2 = 0;
            for (int m = 0; m < DCTCoeffs.GetLength(0); m++)
                for (int b = 0; b < DCTCoeffs.GetLength(1); b++)
                    for (int c = 1; c < DCTCoeffs.GetLength(2); c++)
                        if (Math.Abs(DCTCoeffs[m, b, c]) == 1)
                            counter1++;
                        else if (Math.Abs(DCTCoeffs[m, b, c]) > 1 && Math.Abs(DCTCoeffs[m, b, c]) < 8)
                            counter2++;
            int allBlockCount = Math.Min(counter1 / 15, counter2 / 7);

            int[,] coeffs = SelectDCTCoeffs(DCTCoeffs, allBlockCount);
            int[,] allBlocks = CollectContainerBlocks(coeffs, allBlockCount);
            int[,] filledBlocks = SearchFilledBlocks(allBlocks, length);
            int copiesCount = filledBlocks.GetLength(1) / blockCount;
            (byte[], int[]) messageCopies = ExtractMessageCopies(filledBlocks, length, copiesCount);
            byte[] messageBytes = DecodeMessage(messageCopies.Item1, messageCopies.Item2, length, copiesCount);
            string message = Encoding.ASCII.GetString(messageBytes);
            return message;
        }
        /// <summary>
        /// Встраивание сообщения в изображение
        /// </summary>
        /// <param name="message">Строка сообщения</param>
        /// <param name="newFilename">Имя нового файла, содержащего сообщение</param>
        /// <param name="copiesCount">Количество копий сообщения</param>
        public void EmbedMessage(string message, string newFilename, int copiesCount)
        {
            short[,,] DCTCoeffs = image.DCTCoeffs;
            int counter1 = 0;
            int counter2 = 0;
            int blockCount = 0;
            for (int m = 0; m < DCTCoeffs.GetLength(0); m++)
                for (int b = 0; b < DCTCoeffs.GetLength(1); b++)
                    for (int c = 1; c < DCTCoeffs.GetLength(2); c++)
                        if (Math.Abs(DCTCoeffs[m, b, c]) == 1)
                            counter1++;
                        else if (Math.Abs(DCTCoeffs[m, b, c]) > 1 && Math.Abs(DCTCoeffs[m, b, c]) < 8)
                            counter2++;
            // Количество доступных блоков рабочей области
            maxBlockCount = Math.Min(counter1 / 15, counter2 / 7);
            // Размер сообщения в блоках
            blockCount = (int)Math.Ceiling(message.Length * 8.0 / 11.0) + 4;
            // Извлекаем пригодные для встраивания коэффициенты
            int[,] coeffs = SelectDCTCoeffs(DCTCoeffs, maxBlockCount);
            // Формируем из коэффициентов бинарные блоки контейнера
            int[,] containerBlocks = CollectContainerBlocks(coeffs, blockCount, copiesCount);
            byte[] messageBytes = MessageBytes(message);
            // Формируем бинарные блоки сообщения
            int[,] messageBlocks = GetMessageBlocks(messageBytes, blockCount, copiesCount);
            // Встраиваем синдромы сообщения
            int[,] filledBlocks = EmbedSyndromes(containerBlocks, messageBlocks);
            // Формируем новые коэффициенты из битов стего-контейнера
            int[,] newCoeffs = CoeffsFromContainerBlocks(filledBlocks, coeffs, blockCount, copiesCount);
            // Помещаем новые коэффициенты на свои места
            short[,,] newDCTCoeffs = ReplaceDCTCoeffs(DCTCoeffs, newCoeffs);
            // Записываем новые данные в изображение
            image.WriteNewData(newDCTCoeffs, newFilename);
        }
    }    
}
