﻿using MzLibUtil;

namespace Readers
{
    public enum SupportedFileType
    {
        Ms1Feature,
        Ms2Feature,
        Mzrt_TopFd,
        Ms1Tsv_FlashDeconv,
        Tsv_FlashDeconv,
        ThermoRaw,
        MzML,
        Mgf,
        BrukerD
    }

    public static class SupportedFileTypeExtensions
    {
        /// <summary>
        /// Returns the extension for the file type
        /// </summary>
        /// <param name="type"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        public static string GetFileExtension(this SupportedFileType type)
        {
            return type switch
            {
                SupportedFileType.Ms1Feature => "_ms1.feature",
                SupportedFileType.Ms2Feature => "_ms2.feature",
                SupportedFileType.Mzrt_TopFd => ".mzrt.csv",
                SupportedFileType.Ms1Tsv_FlashDeconv => "_ms1.tsv",
                SupportedFileType.Tsv_FlashDeconv => ".tsv",
                SupportedFileType.ThermoRaw => ".raw",
                SupportedFileType.MzML => ".mzML",
                SupportedFileType.Mgf => ".mgf",
                SupportedFileType.BrukerD => ".d",
                _ => throw new MzLibException("File type not supported")
            };
        }

        public static SupportedFileType ParseFileType(this string filePath)
        {
            switch (Path.GetExtension(filePath).ToLower())
            {
                case ".raw": return SupportedFileType.ThermoRaw;
                case ".mzml": return SupportedFileType.MzML;
                case ".mgf": return SupportedFileType.Mgf;
                case ".d": return SupportedFileType.BrukerD;
                case ".feature":
                    if (filePath.EndsWith(SupportedFileType.Ms1Feature.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Ms1Feature;
                    if (filePath.EndsWith(SupportedFileType.Ms2Feature.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Ms2Feature;
                    throw new MzLibException("Feature file type not supported");

                case ".csv":
                    if (filePath.EndsWith(SupportedFileType.Mzrt_TopFd.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Mzrt_TopFd;
                    throw new MzLibException("Csv file type not supported");

                case ".tsv":
                    if (filePath.EndsWith(SupportedFileType.Ms1Tsv_FlashDeconv.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Ms1Tsv_FlashDeconv;
                    // catchall for other tsv types, one one implemented right now
                    if (filePath.EndsWith(SupportedFileType.Tsv_FlashDeconv.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase) &&
                        !filePath.EndsWith(SupportedFileType.Ms1Tsv_FlashDeconv.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Tsv_FlashDeconv;
                    throw new MzLibException("Tsv file type not supported");

                default:
                    throw new MzLibException("File type not supported");
            }
        }
    }
}
