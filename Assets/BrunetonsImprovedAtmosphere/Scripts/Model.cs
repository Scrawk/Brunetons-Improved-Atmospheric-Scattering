using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;

namespace BrunetonsImprovedAtmosphere
{

    public class Model
    {
        private const int READ = 0;
        private const int WRITE = 1;

        private const double kLambdaR = 680.0;
        private const double kLambdaG = 550.0;
        private const double kLambdaB = 440.0;

        private const int kLambdaMin = 360;
        private const int kLambdaMax = 830;

        /// <summary>
        /// The wavelength values, in nanometers, and sorted in increasing order, for
        /// which the solar_irradiance, rayleigh_scattering, mie_scattering,
        /// mie_extinction and ground_albedo samples are provided. If your shaders
        /// use luminance values (as opposed to radiance values, see above), use a
        /// large number of wavelengths (e.g. between 15 and 50) to get accurate
        /// results (this number of wavelengths has absolutely no impact on the
        /// shader performance).
        /// </summary>
        public IList<double> Wavelengths { get; set; }

        /// <summary>
        /// The solar irradiance at the top of the atmosphere, in W/m^2/nm. This
        /// vector must have the same size as the wavelengths parameter.
        /// </summary>
        public IList<double> SolarIrradiance { get; set; }

        /// <summary>
        /// The sun's angular radius, in radians. Warning: the implementation uses
        /// approximations that are valid only if this value is smaller than 0.1.
        /// </summary>
        public double SunAngularRadius { get; set; }

        /// <summary>
        /// The distance between the planet center and the bottom of the atmosphere in m.
        /// </summary>
        public double BottomRadius { get; set; }

        /// <summary>
        /// The distance between the planet center and the top of the atmosphere in m.
        /// </summary>
        public double TopRadius { get; set; }

        /// <summary>
        /// The density profile of air molecules, i.e. a function from altitude to
        /// dimensionless values between 0 (null density) and 1 (maximum density).
        /// Layers must be sorted from bottom to top. The width of the last layer is
        /// ignored, i.e. it always extend to the top atmosphere boundary. At most 2
        /// layers can be specified.
        /// </summary>
        public DensityProfileLayer RayleighDensity { get; set; }

        /// <summary>
        /// The scattering coefficient of air molecules at the altitude where their
        /// density is maximum (usually the bottom of the atmosphere), as a function
        /// of wavelength, in m^-1. The scattering coefficient at altitude h is equal
        /// to 'rayleigh_scattering' times 'rayleigh_density' at this altitude. This
        /// vector must have the same size as the wavelengths parameter.
        /// </summary>
        public IList<double> RayleighScattering { get; set; }

        /// <summary>
        /// The density profile of aerosols, i.e. a function from altitude to
        /// dimensionless values between 0 (null density) and 1 (maximum density).
        /// Layers must be sorted from bottom to top. The width of the last layer is
        /// ignored, i.e. it always extend to the top atmosphere boundary. At most 2
        /// layers can be specified.
        /// </summary>
        public DensityProfileLayer MieDensity { get; set; }

        /// <summary>
        /// The scattering coefficient of aerosols at the altitude where their
        /// density is maximum (usually the bottom of the atmosphere), as a function
        /// of wavelength, in m^-1. The scattering coefficient at altitude h is equal
        /// to 'mie_scattering' times 'mie_density' at this altitude. This vector
        /// must have the same size as the wavelengths parameter.
        /// </summary>
        public IList<double> MieScattering { get; set; }

        /// <summary>
        /// The extinction coefficient of aerosols at the altitude where their
        /// density is maximum (usually the bottom of the atmosphere), as a function
        /// of wavelength, in m^-1. The extinction coefficient at altitude h is equal
        /// to 'mie_extinction' times 'mie_density' at this altitude. This vector
        /// must have the same size as the wavelengths parameter.
        /// </summary>
        public IList<double> MieExtinction { get; set; }

        /// <summary>
        /// The asymetry parameter for the Cornette-Shanks phase function for the aerosols.
        /// </summary>
        public double MiePhaseFunctionG { get; set; }

        /// <summary>
        /// The density profile of air molecules that absorb light (e.g. ozone), i.e.
        /// a function from altitude to dimensionless values between 0 (null density)
        /// and 1 (maximum density). Layers must be sorted from bottom to top. The
        /// width of the last layer is ignored, i.e. it always extend to the top
        /// atmosphere boundary. At most 2 layers can be specified.
        /// </summary>
        public IList<DensityProfileLayer> AbsorptionDensity { get; set; }

        /// <summary>
        /// The extinction coefficient of molecules that absorb light (e.g. ozone) at
        /// the altitude where their density is maximum, as a function of wavelength,
        /// in m^-1. The extinction coefficient at altitude h is equal to
        /// 'absorption_extinction' times 'absorption_density' at this altitude. This
        /// vector must have the same size as the wavelengths parameter.
        /// </summary>
        public IList<double> AbsorptionExtinction { get; set; }

        /// <summary>
        /// The average albedo of the ground, as a function of wavelength. This
        /// vector must have the same size as the wavelengths parameter.
        /// </summary>
        public IList<double> GroundAlbedo { get; set; }

        /// <summary>
        /// The maximum Sun zenith angle for which atmospheric scattering must be
        /// precomputed, in radians (for maximum precision, use the smallest Sun
        /// zenith angle yielding negligible sky light radiance values. For instance,
        /// for the Earth case, 102 degrees is a good choice for most cases (120
        /// degrees is necessary for very high exposure values).
        /// </summary>
        public double MaxSunZenithAngle { get; set; }

        /// <summary>
        /// The length unit used in your shaders and meshes. This is the length unit
        /// which must be used when calling the atmosphere model shader functions.
        /// </summary>
        public double LengthUnitInMeters { get; set; }

        /// <summary>
        /// The number of wavelengths for which atmospheric scattering must be
        /// precomputed (the temporary GPU memory used during precomputations, and
        /// the GPU memory used by the precomputed results, is independent of this
        /// number, but the precomputation time is directly proportional to this number):
        /// - if this number is less than or equal to 3, scattering is precomputed
        /// for 3 wavelengths, and stored as irradiance values. Then both the
        /// radiance-based and the luminance-based API functions are provided (see
        /// the above note).
        /// - otherwise, scattering is precomputed for this number of wavelengths
        /// (rounded up to a multiple of 3), integrated with the CIE color matching
        /// functions, and stored as illuminance values. Then only the
        /// luminance-based API functions are provided (see the above note).
        /// </summary>
        public int NumPrecomputedWavelengths { get { return UseLuminance == LUMINANCE.PRECOMPUTED ? 15 : 3; } }

        /// <summary>
        /// Whether to pack the (red component of the) single Mie scattering with the
        /// Rayleigh and multiple scattering in a single texture, or to store the
        /// (3 components of the) single Mie scattering in a separate texture.
        /// </summary>
        public bool CombineScatteringTextures { get; set; }

        /// <summary>
        /// Use radiance or luminance mode.
        /// </summary>
        public LUMINANCE UseLuminance { get; set; }

        /// <summary>
        /// Whether to use half precision floats (16 bits) or single precision floats
        /// (32 bits) for the precomputed textures. Half precision is sufficient for
        /// most cases, except for very high exposure values.
        /// </summary>
        public bool HalfPrecision { get; set; }

        public RenderTexture TransmittanceTexture { get; private set; }

        public RenderTexture ScatteringTexture { get; private set; }

        public RenderTexture IrradianceTexture { get; private set; }

        public RenderTexture OptionalSingleMieScatteringTexture { get; private set; }

        public Model()
        {

        }

        /// <summary>
        /// Bind to a pixel shader for rendering.
        /// </summary>
        public void BindToMaterial(Material mat)
        {

            if (UseLuminance == LUMINANCE.NONE)
                mat.EnableKeyword("RADIANCE_API_ENABLED");
            else
                mat.DisableKeyword("RADIANCE_API_ENABLED");

            if (CombineScatteringTextures)
                mat.EnableKeyword("COMBINED_SCATTERING_TEXTURES");
            else
                mat.DisableKeyword("COMBINED_SCATTERING_TEXTURES");

            mat.SetTexture("transmittance_texture", TransmittanceTexture);
            mat.SetTexture("scattering_texture", ScatteringTexture);
            mat.SetTexture("irradiance_texture", IrradianceTexture);

            if(CombineScatteringTextures)
                mat.SetTexture("single_mie_scattering_texture", Texture2D.blackTexture);
            else
                mat.SetTexture("single_mie_scattering_texture", OptionalSingleMieScatteringTexture);

            mat.SetInt("TRANSMITTANCE_TEXTURE_WIDTH", CONSTANTS.TRANSMITTANCE_WIDTH);
            mat.SetInt("TRANSMITTANCE_TEXTURE_HEIGHT", CONSTANTS.TRANSMITTANCE_HEIGHT);
            mat.SetInt("SCATTERING_TEXTURE_R_SIZE", CONSTANTS.SCATTERING_R);
            mat.SetInt("SCATTERING_TEXTURE_MU_SIZE", CONSTANTS.SCATTERING_MU);
            mat.SetInt("SCATTERING_TEXTURE_MU_S_SIZE", CONSTANTS.SCATTERING_MU_S);
            mat.SetInt("SCATTERING_TEXTURE_NU_SIZE", CONSTANTS.SCATTERING_NU);
            mat.SetInt("SCATTERING_TEXTURE_WIDTH", CONSTANTS.SCATTERING_WIDTH);
            mat.SetInt("SCATTERING_TEXTURE_HEIGHT", CONSTANTS.SCATTERING_HEIGHT);
            mat.SetInt("SCATTERING_TEXTURE_DEPTH", CONSTANTS.SCATTERING_DEPTH);
            mat.SetInt("IRRADIANCE_TEXTURE_WIDTH", CONSTANTS.IRRADIANCE_WIDTH);
            mat.SetInt("IRRADIANCE_TEXTURE_HEIGHT", CONSTANTS.IRRADIANCE_HEIGHT);

            mat.SetFloat("sun_angular_radius", (float)SunAngularRadius);
            mat.SetFloat("bottom_radius", (float)(BottomRadius / LengthUnitInMeters));
            mat.SetFloat("top_radius", (float)(TopRadius / LengthUnitInMeters));
            mat.SetFloat("mie_phase_function_g", (float)MiePhaseFunctionG);
            mat.SetFloat("mu_s_min", (float)Math.Cos(MaxSunZenithAngle));

            Vector3 skySpectralRadianceToLuminance, sunSpectralRadianceToLuminance;
            SkySunRadianceToLuminance(out skySpectralRadianceToLuminance, out sunSpectralRadianceToLuminance);

            mat.SetVector("SKY_SPECTRAL_RADIANCE_TO_LUMINANCE", skySpectralRadianceToLuminance);
            mat.SetVector("SUN_SPECTRAL_RADIANCE_TO_LUMINANCE", sunSpectralRadianceToLuminance);

            double[] lambdas = new double[] { kLambdaR, kLambdaG, kLambdaB };

            Vector3 solarIrradiance = ToVector(Wavelengths, SolarIrradiance, lambdas, 1.0);
            mat.SetVector("solar_irradiance", solarIrradiance);

            Vector3 rayleighScattering = ToVector(Wavelengths, RayleighScattering, lambdas, LengthUnitInMeters);
            mat.SetVector("rayleigh_scattering", rayleighScattering);

            Vector3 mieScattering = ToVector(Wavelengths, MieScattering, lambdas, LengthUnitInMeters);
            mat.SetVector("mie_scattering", mieScattering);
        }

        public void Release()
        {
            ReleaseTexture(TransmittanceTexture);
            ReleaseTexture(ScatteringTexture);
            ReleaseTexture(IrradianceTexture);
            ReleaseTexture(OptionalSingleMieScatteringTexture);
        }

        /// <summary>
        /// The Init method precomputes the atmosphere textures. It first allocates the
        /// temporary resources it needs, then calls Precompute to do the
        /// actual precomputations, and finally destroys the temporary resources.
        ///
        /// Note that there are two precomputation modes here, depending on whether we
        /// want to store precomputed irradiance or illuminance values:
        ///
        /// In precomputed irradiance mode, we simply need to call
        /// Precompute with the 3 wavelengths for which we want to precompute
        /// irradiance, namely kLambdaR, kLambdaG, kLambdaB(with the identity matrix for
        /// luminance_from_radiance, since we don't want any conversion from radiance to luminance).
        /// 
        /// In precomputed illuminance mode, we need to precompute irradiance for
        /// num_precomputed_wavelengths, and then integrate the results,
        /// multiplied with the 3 CIE xyz color matching functions and the XYZ to sRGB
        /// matrix to get sRGB illuminance values.
        /// A naive solution would be to allocate temporary textures for the
        /// intermediate irradiance results, then perform the integration from irradiance
        /// to illuminance and store the result in the final precomputed texture. In
        /// pseudo-code (and assuming one wavelength per texture instead of 3):
        ///  
        ///  create n temporary irradiance textures
        ///  for each wavelength lambda in the n wavelengths:
        ///     precompute irradiance at lambda into one of the temporary textures
        ///  initializes the final illuminance texture with zeros
        ///  for each wavelength lambda in the n wavelengths:
        ///     accumulate in the final illuminance texture the product of the
        ///     precomputed irradiance at lambda (read from the temporary textures)
        ///     with the value of the 3 sRGB color matching functions at lambda 
        ///     (i.e. the product of the XYZ to sRGB matrix with the CIE xyz color matching functions).
        ///  
        /// However, this be would waste GPU memory. Instead, we can avoid allocating
        /// temporary irradiance textures, by merging the two above loops:
        ///  
        ///   for each wavelength lambda in the n wavelengths:
        ///     accumulate in the final illuminance texture (or, for the first
        ///     iteration, set this texture to) the product of the precomputed
        ///     irradiance at lambda (computed on the fly) with the value of the 3
        ///     sRGB color matching functions at lambda.
        ///  
        /// This is the method we use below, with 3 wavelengths per iteration instead
        /// of 1, using Precompute to compute 3 irradiances values per
        /// iteration, and luminance_from_radiance to multiply 3 irradiances
        /// with the values of the 3 sRGB color matching functions at 3 different
        /// wavelengths (yielding a 3x3 matrix).
        ///
        /// This yields the following implementation:
        /// </summary>
        public void Init(ComputeShader compute, int num_scattering_orders)
        {

            // The precomputations require temporary textures, in particular to store the
            // contribution of one scattering order, which is needed to compute the next
            // order of scattering (the final precomputed textures store the sum of all
            // the scattering orders). We allocate them here, and destroy them at the end
            // of this method.
            TextureBuffer buffer = new TextureBuffer(HalfPrecision);
            buffer.Clear(compute);

            // The actual precomputations depend on whether we want to store precomputed
            // irradiance or illuminance values.
            if (NumPrecomputedWavelengths <= 3) 
            {
                Precompute(compute, buffer, null, null, false, num_scattering_orders);
            } 
            else 
            {
                int num_iterations = (NumPrecomputedWavelengths + 2) / 3;
                double dlambda = (kLambdaMax - kLambdaMin) / (3.0 * num_iterations);

                for (int i = 0; i < num_iterations; ++i)
                {
                    double[] lambdas = new double[]
                    {
                        kLambdaMin + (3 * i + 0.5) * dlambda,
                        kLambdaMin + (3 * i + 1.5) * dlambda,
                        kLambdaMin + (3 * i + 2.5) * dlambda
                    };

                    double[] luminance_from_radiance = new double[]
                    {
                        Coeff(lambdas[0], 0) * dlambda, Coeff(lambdas[1], 0) * dlambda, Coeff(lambdas[2], 0) * dlambda,
                        Coeff(lambdas[0], 1) * dlambda, Coeff(lambdas[1], 1) * dlambda, Coeff(lambdas[2], 1) * dlambda,
                        Coeff(lambdas[0], 2) * dlambda, Coeff(lambdas[1], 2) * dlambda, Coeff(lambdas[2], 2) * dlambda
                    };

                    bool blend = i > 0;
                    Precompute(compute, buffer, lambdas, luminance_from_radiance, blend, num_scattering_orders);
                }

                // After the above iterations, the transmittance texture contains the
                // transmittance for the 3 wavelengths used at the last iteration. But we
                // want the transmittance at kLambdaR, kLambdaG, kLambdaB instead, so we
                // must recompute it here for these 3 wavelengths:
                int compute_transmittance = compute.FindKernel("ComputeTransmittance");
                BindToCompute(compute, null, null);
                compute.SetTexture(compute_transmittance, "transmittanceWrite", buffer.TransmittanceArray[WRITE]);
                compute.SetVector("blend", new Vector4(0, 0, 0, 0));

                int NUM = CONSTANTS.NUM_THREADS;
                compute.Dispatch(compute_transmittance, CONSTANTS.TRANSMITTANCE_WIDTH / NUM, CONSTANTS.TRANSMITTANCE_HEIGHT / NUM, 1);
                Swap(buffer.TransmittanceArray);
            }

            //Grab ref to textures and mark as null in buffer so they are not released.
            TransmittanceTexture = buffer.TransmittanceArray[READ];
            buffer.TransmittanceArray[READ] = null;

            ScatteringTexture = buffer.ScatteringArray[READ];
            buffer.ScatteringArray[READ] = null;

            IrradianceTexture = buffer.IrradianceArray[READ];
            buffer.IrradianceArray[READ] = null;

            if(CombineScatteringTextures)
            {
                OptionalSingleMieScatteringTexture = null;
            }
            else
            {
                OptionalSingleMieScatteringTexture = buffer.OptionalSingleMieScatteringArray[READ];
                buffer.OptionalSingleMieScatteringArray[READ] = null;
            }

            // Delete the temporary resources allocated at the begining of this method.
            buffer.Release();

        }

        private double Coeff(double lambda, int component) 
        {
            // Note that we don't include MAX_LUMINOUS_EFFICACY here, to avoid
            // artefacts due to too large values when using half precision on GPU.
            // We add this term back in kAtmosphereShader, via
            // SKY_SPECTRAL_RADIANCE_TO_LUMINANCE (see also the comments in the
            // Model constructor).
            double x = CieColorMatchingFunctionTableValue(lambda, 1);
            double y = CieColorMatchingFunctionTableValue(lambda, 2);
            double z = CieColorMatchingFunctionTableValue(lambda, 3);
            double sRGB = CONSTANTS.XYZ_TO_SRGB[component * 3 + 0] * x + 
                          CONSTANTS.XYZ_TO_SRGB[component * 3 + 1] * y + 
                          CONSTANTS.XYZ_TO_SRGB[component * 3 + 2] * z;

            return sRGB;
        }

        /// <summary>
        /// Bind to a compute shader for precomutation of textures.
        /// </summary>
        private void BindToCompute(ComputeShader compute, double[] lambdas, double[] luminance_from_radiance)
        {
            if(lambdas == null)
                lambdas = new double[] { kLambdaR, kLambdaG, kLambdaB };

            if(luminance_from_radiance == null)
                luminance_from_radiance = new double[] { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

            compute.SetInt("TRANSMITTANCE_TEXTURE_WIDTH", CONSTANTS.TRANSMITTANCE_WIDTH);
            compute.SetInt("TRANSMITTANCE_TEXTURE_HEIGHT", CONSTANTS.TRANSMITTANCE_HEIGHT);
            compute.SetInt("SCATTERING_TEXTURE_R_SIZE", CONSTANTS.SCATTERING_R);
            compute.SetInt("SCATTERING_TEXTURE_MU_SIZE", CONSTANTS.SCATTERING_MU);
            compute.SetInt("SCATTERING_TEXTURE_MU_S_SIZE", CONSTANTS.SCATTERING_MU_S);
            compute.SetInt("SCATTERING_TEXTURE_NU_SIZE", CONSTANTS.SCATTERING_NU);
            compute.SetInt("SCATTERING_TEXTURE_WIDTH", CONSTANTS.SCATTERING_WIDTH);
            compute.SetInt("SCATTERING_TEXTURE_HEIGHT", CONSTANTS.SCATTERING_HEIGHT);
            compute.SetInt("SCATTERING_TEXTURE_DEPTH", CONSTANTS.SCATTERING_DEPTH);
            compute.SetInt("IRRADIANCE_TEXTURE_WIDTH", CONSTANTS.IRRADIANCE_WIDTH);
            compute.SetInt("IRRADIANCE_TEXTURE_HEIGHT", CONSTANTS.IRRADIANCE_HEIGHT);

            Vector3 skySpectralRadianceToLuminance, sunSpectralRadianceToLuminance;
            SkySunRadianceToLuminance(out skySpectralRadianceToLuminance, out sunSpectralRadianceToLuminance);

            compute.SetVector("SKY_SPECTRAL_RADIANCE_TO_LUMINANCE", skySpectralRadianceToLuminance);
            compute.SetVector("SUN_SPECTRAL_RADIANCE_TO_LUMINANCE", sunSpectralRadianceToLuminance);

            Vector3 solarIrradiance = ToVector(Wavelengths, SolarIrradiance, lambdas, 1.0);
            compute.SetVector("solar_irradiance", solarIrradiance);

            Vector3 rayleighScattering = ToVector(Wavelengths, RayleighScattering, lambdas, LengthUnitInMeters);
            BindDensityLayer(compute, RayleighDensity);
            compute.SetVector("rayleigh_scattering", rayleighScattering);

            Vector3 mieScattering = ToVector(Wavelengths, MieScattering, lambdas, LengthUnitInMeters);
            Vector3 mieExtinction = ToVector(Wavelengths, MieExtinction, lambdas, LengthUnitInMeters);
            BindDensityLayer(compute, MieDensity);
            compute.SetVector("mie_scattering", mieScattering);
            compute.SetVector("mie_extinction", mieExtinction);

            Vector3 absorptionExtinction = ToVector(Wavelengths, AbsorptionExtinction, lambdas, LengthUnitInMeters);
            BindDensityLayer(compute, AbsorptionDensity[0]);
            BindDensityLayer(compute, AbsorptionDensity[1]);
            compute.SetVector("absorption_extinction", absorptionExtinction);

            Vector3 groundAlbedo = ToVector(Wavelengths, GroundAlbedo, lambdas, 1.0);
            compute.SetVector("ground_albedo", groundAlbedo);

            compute.SetFloats("luminanceFromRadiance", ToMatrix(luminance_from_radiance));
            compute.SetFloat("sun_angular_radius", (float)SunAngularRadius);
            compute.SetFloat("bottom_radius", (float)(BottomRadius / LengthUnitInMeters));
            compute.SetFloat("top_radius", (float)(TopRadius / LengthUnitInMeters));
            compute.SetFloat("mie_phase_function_g", (float)MiePhaseFunctionG);
            compute.SetFloat("mu_s_min", (float)Math.Cos(MaxSunZenithAngle));
        }

        private void BindDensityLayer(ComputeShader compute, DensityProfileLayer layer)
        {
            compute.SetFloat(layer.Name + "_width", (float)(layer.Width / LengthUnitInMeters));
            compute.SetFloat(layer.Name + "_exp_term", (float)layer.ExpTerm);
            compute.SetFloat(layer.Name + "_exp_scale", (float)(layer.ExpScale * LengthUnitInMeters));
            compute.SetFloat(layer.Name + "_linear_term", (float)(layer.LinearTerm * LengthUnitInMeters));
            compute.SetFloat(layer.Name + "_constant_term", (float)layer.ConstantTerm);
        }

        private Vector3 ToVector(IList<double> wavelengths, IList<double> v, IList<double> lambdas, double scale)
        {
            double r = Interpolate(wavelengths, v, lambdas[0]) * scale;
            double g = Interpolate(wavelengths, v, lambdas[1]) * scale;
            double b = Interpolate(wavelengths, v, lambdas[2]) * scale;

            return new Vector3((float)r, (float)g, (float)b);
        }

        /// <summary>
        /// Finally, we need a utility function to compute the value of the conversion
        /// constants *_RADIANCE_TO_LUMINANCE, used above to convert the
        /// spectral results into luminance values. These are the constants k_r, k_g, k_b
        /// described in Section 14.3 of <a href="https://arxiv.org/pdf/1612.04336.pdf">A
        /// Qualitative and Quantitative Evaluation of 8 Clear Sky Models</a>.
        ///
        /// Computing their value requires an integral of a function times a CIE color
        /// matching function. Thus, we first need functions to interpolate an arbitrary
        /// function (specified by some samples), and a CIE color matching function
        /// (specified by tabulated values), at an arbitrary wavelength. This is the purpose
        /// of the following two functions:
        /// </summary>
        private static double CieColorMatchingFunctionTableValue(double wavelength, int column)
        {
            if (wavelength <= kLambdaMin || wavelength >= kLambdaMax) return 0.0;

            double u = (wavelength - kLambdaMin) / 5.0;
            int row = (int)Math.Floor(u);

            u -= row;
            return CONSTANTS.CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[4 * row + column] * (1.0 - u) + CONSTANTS.CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[4 * (row + 1) + column] * u;
        }

        private static double Interpolate(IList<double> wavelengths, IList<double> wavelength_function, double wavelength)
        {
            if (wavelength < wavelengths[0]) return wavelength_function[0];

            for (int i = 0; i < wavelengths.Count - 1; ++i)
            {
                if (wavelength < wavelengths[i + 1])
                {
                    double u = (wavelength - wavelengths[i]) / (wavelengths[i + 1] - wavelengths[i]);
                    return wavelength_function[i] * (1.0 - u) + wavelength_function[i + 1] * u;
                }
            }

            return wavelength_function[wavelength_function.Count - 1];
        }

        /// <summary>
        ///  Compute the values for the SKY_RADIANCE_TO_LUMINANCE constant. In theory
        /// this should be 1 in precomputed illuminance mode (because the precomputed
        /// textures already contain illuminance values). In practice, however, storing
        /// true illuminance values in half precision textures yields artefacts
        /// (because the values are too large), so we store illuminance values divided
        /// by MAX_LUMINOUS_EFFICACY instead. This is why, in precomputed illuminance
        /// mode, we set SKY_RADIANCE_TO_LUMINANCE to MAX_LUMINOUS_EFFICACY.
        /// </summary>
        private void SkySunRadianceToLuminance(out Vector3 skySpectralRadianceToLuminance, out Vector3 sunSpectralRadianceToLuminance)
        {
            bool precompute_illuminance = NumPrecomputedWavelengths > 3;
            double sky_k_r, sky_k_g, sky_k_b;

            if (precompute_illuminance)
                sky_k_r = sky_k_g = sky_k_b = CONSTANTS.MAX_LUMINOUS_EFFICACY;
            else
                ComputeSpectralRadianceToLuminanceFactors(Wavelengths, SolarIrradiance, -3, out sky_k_r, out sky_k_g, out sky_k_b);

            // Compute the values for the SUN_RADIANCE_TO_LUMINANCE constant.
            double sun_k_r, sun_k_g, sun_k_b;
            ComputeSpectralRadianceToLuminanceFactors(Wavelengths, SolarIrradiance, 0, out sun_k_r, out sun_k_g, out sun_k_b);

            skySpectralRadianceToLuminance = new Vector3((float)sky_k_r, (float)sky_k_g, (float)sky_k_b);
            sunSpectralRadianceToLuminance = new Vector3((float)sun_k_r, (float)sun_k_g, (float)sun_k_b);
        }

        /// <summary>
        /// We can then implement a utility function to compute the "spectral radiance to
        /// luminance" conversion constants (see Section 14.3 in <a
        /// href="https://arxiv.org/pdf/1612.04336.pdf">A Qualitative and Quantitative
        /// Evaluation of 8 Clear Sky Models</a> for their definitions):
        /// The returned constants are in lumen.nm / watt.
        /// </summary>
        private static void ComputeSpectralRadianceToLuminanceFactors(IList<double> wavelengths, IList<double> solar_irradiance, double lambda_power, out double k_r, out double k_g, out double k_b) 
        {
            k_r = 0.0;
            k_g = 0.0;
            k_b = 0.0;
            double solar_r = Interpolate(wavelengths, solar_irradiance, kLambdaR);
            double solar_g = Interpolate(wavelengths, solar_irradiance, kLambdaG);
            double solar_b = Interpolate(wavelengths, solar_irradiance, kLambdaB);
            int dlambda = 1;

            for (int lambda = kLambdaMin; lambda < kLambdaMax; lambda += dlambda) 
            {
                double x_bar = CieColorMatchingFunctionTableValue(lambda, 1);
                double y_bar = CieColorMatchingFunctionTableValue(lambda, 2);
                double z_bar = CieColorMatchingFunctionTableValue(lambda, 3);

                double[] xyz2srgb = CONSTANTS.XYZ_TO_SRGB;
                double r_bar = xyz2srgb[0] * x_bar + xyz2srgb[1] * y_bar + xyz2srgb[2] * z_bar;
                double g_bar = xyz2srgb[3] * x_bar + xyz2srgb[4] * y_bar + xyz2srgb[5] * z_bar;
                double b_bar = xyz2srgb[6] * x_bar + xyz2srgb[7] * y_bar + xyz2srgb[8] * z_bar;
                double irradiance = Interpolate(wavelengths, solar_irradiance, lambda);

                k_r += r_bar * irradiance / solar_r * Math.Pow(lambda / kLambdaR, lambda_power);
                k_g += g_bar * irradiance / solar_g * Math.Pow(lambda / kLambdaG, lambda_power);
                k_b += b_bar * irradiance / solar_b * Math.Pow(lambda / kLambdaB, lambda_power);
            }

            k_r *= CONSTANTS.MAX_LUMINOUS_EFFICACY * dlambda;
            k_g *= CONSTANTS.MAX_LUMINOUS_EFFICACY * dlambda;
            k_b *= CONSTANTS.MAX_LUMINOUS_EFFICACY * dlambda;
        }

        /// <summary>
        /// The utility method ConvertSpectrumToLinearSrgb is implemented
        /// with a simple numerical integration of the given function, times the CIE color
        /// matching funtions(with an integration step of 1nm), followed by a matrix
        /// multiplication:
        /// </summary>
        public void ConvertSpectrumToLinearSrgb(out double r, out double g, out double b)
        {
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            const int dlambda = 1;
            for (int lambda = kLambdaMin; lambda < kLambdaMax; lambda += dlambda)
            {
                double value = Interpolate(Wavelengths, SolarIrradiance, lambda);
                x += CieColorMatchingFunctionTableValue(lambda, 1) * value;
                y += CieColorMatchingFunctionTableValue(lambda, 2) * value;
                z += CieColorMatchingFunctionTableValue(lambda, 3) * value;
            }

            double[] XYZ_TO_SRGB = CONSTANTS.XYZ_TO_SRGB;
            r = CONSTANTS.MAX_LUMINOUS_EFFICACY * (XYZ_TO_SRGB[0] * x + XYZ_TO_SRGB[1] * y + XYZ_TO_SRGB[2] * z) * dlambda;
            g = CONSTANTS.MAX_LUMINOUS_EFFICACY * (XYZ_TO_SRGB[3] * x + XYZ_TO_SRGB[4] * y + XYZ_TO_SRGB[5] * z) * dlambda;
            b = CONSTANTS.MAX_LUMINOUS_EFFICACY * (XYZ_TO_SRGB[6] * x + XYZ_TO_SRGB[7] * y + XYZ_TO_SRGB[8] * z) * dlambda;
        }

        /// <summary>
        /// Finally, we provide the actual implementation of the precomputation algorithm
        /// described in Algorithm 4.1 of
        /// <a href="https://hal.inria.fr/inria-00288758/en">our paper</a>. Each step is
        /// explained by the inline comments below.
        /// </summary>
        void Precompute(
            ComputeShader compute,
            TextureBuffer buffer,
            double[] lambdas,
            double[] luminance_from_radiance,
            bool blend,
            int num_scattering_orders)
        {

            int BLEND = blend ? 1 : 0;
            int NUM_THREADS = CONSTANTS.NUM_THREADS;

            BindToCompute(compute, lambdas, luminance_from_radiance);

            int compute_transmittance = compute.FindKernel("ComputeTransmittance");
            int compute_direct_irradiance = compute.FindKernel("ComputeDirectIrradiance");
            int compute_single_scattering = compute.FindKernel("ComputeSingleScattering");
            int compute_scattering_density = compute.FindKernel("ComputeScatteringDensity");
            int compute_indirect_irradiance = compute.FindKernel("ComputeIndirectIrradiance");
            int compute_multiple_scattering = compute.FindKernel("ComputeMultipleScattering");
 
            // Compute the transmittance, and store it in transmittance_texture
            compute.SetTexture(compute_transmittance, "transmittanceWrite", buffer.TransmittanceArray[WRITE]);
            compute.SetVector("blend", new Vector4(0, 0, 0, 0));
            compute.Dispatch(compute_transmittance, CONSTANTS.TRANSMITTANCE_WIDTH / NUM_THREADS, CONSTANTS.TRANSMITTANCE_HEIGHT / NUM_THREADS, 1);
            Swap(buffer.TransmittanceArray);

            // Compute the direct irradiance, store it in delta_irradiance_texture and,
            // depending on 'blend', either initialize irradiance_texture_ with zeros or
            // leave it unchanged (we don't want the direct irradiance in
            // irradiance_texture_, but only the irradiance from the sky).
            compute.SetTexture(compute_direct_irradiance, "deltaIrradianceWrite", buffer.DeltaIrradianceTexture); //0
            compute.SetTexture(compute_direct_irradiance, "irradianceWrite", buffer.IrradianceArray[WRITE]); //1
            compute.SetTexture(compute_direct_irradiance, "irradianceRead", buffer.IrradianceArray[READ]);
            compute.SetTexture(compute_direct_irradiance, "transmittanceRead", buffer.TransmittanceArray[READ]);
            compute.SetVector("blend", new Vector4(0, BLEND, 0, 0));
            compute.Dispatch(compute_direct_irradiance, CONSTANTS.IRRADIANCE_WIDTH / NUM_THREADS, CONSTANTS.IRRADIANCE_HEIGHT / NUM_THREADS, 1);
            Swap(buffer.IrradianceArray);

            // Compute the rayleigh and mie single scattering, store them in
            // delta_rayleigh_scattering_texture and delta_mie_scattering_texture, and
            // either store them or accumulate them in scattering_texture_ and
            // optional_single_mie_scattering_texture_.
            compute.SetTexture(compute_single_scattering, "deltaRayleighScatteringWrite", buffer.DeltaRayleighScatteringTexture); //0
            compute.SetTexture(compute_single_scattering, "deltaMieScatteringWrite", buffer.DeltaMieScatteringTexture); //1
            compute.SetTexture(compute_single_scattering, "scatteringWrite", buffer.ScatteringArray[WRITE]); //2
            compute.SetTexture(compute_single_scattering, "scatteringRead", buffer.ScatteringArray[READ]);
            compute.SetTexture(compute_single_scattering, "singleMieScatteringWrite", buffer.OptionalSingleMieScatteringArray[WRITE]); //3
            compute.SetTexture(compute_single_scattering, "singleMieScatteringRead", buffer.OptionalSingleMieScatteringArray[READ]);
            compute.SetTexture(compute_single_scattering, "transmittanceRead", buffer.TransmittanceArray[READ]);
            compute.SetVector("blend", new Vector4(0, 0, BLEND, BLEND));

            for (int layer = 0; layer < CONSTANTS.SCATTERING_DEPTH; ++layer) 
            {
                compute.SetInt("layer", layer);
                compute.Dispatch(compute_single_scattering, CONSTANTS.SCATTERING_WIDTH / NUM_THREADS, CONSTANTS.SCATTERING_HEIGHT / NUM_THREADS, 1);
            }
            Swap(buffer.ScatteringArray);
            Swap(buffer.OptionalSingleMieScatteringArray);

            // Compute the 2nd, 3rd and 4th order of scattering, in sequence.
            for (int scattering_order = 2; scattering_order <= num_scattering_orders; ++scattering_order) 
            {
                // Compute the scattering density, and store it in
                // delta_scattering_density_texture.
                compute.SetTexture(compute_scattering_density, "deltaScatteringDensityWrite", buffer.DeltaScatteringDensityTexture); //0
                compute.SetTexture(compute_scattering_density, "transmittanceRead", buffer.TransmittanceArray[READ]);
                compute.SetTexture(compute_scattering_density, "singleRayleighScatteringRead", buffer.DeltaRayleighScatteringTexture);
                compute.SetTexture(compute_scattering_density, "singleMieScatteringRead", buffer.DeltaMieScatteringTexture);
                compute.SetTexture(compute_scattering_density, "multipleScatteringRead", buffer.DeltaMultipleScatteringTexture);
                compute.SetTexture(compute_scattering_density, "irradianceRead", buffer.DeltaIrradianceTexture);
                compute.SetInt("scatteringOrder", scattering_order);
                compute.SetVector("blend", new Vector4(0, 0, 0, 0));

                for (int layer = 0; layer < CONSTANTS.SCATTERING_DEPTH; ++layer) 
                {
                    compute.SetInt("layer", layer);
                    compute.Dispatch(compute_scattering_density, CONSTANTS.SCATTERING_WIDTH / NUM_THREADS, CONSTANTS.SCATTERING_HEIGHT / NUM_THREADS, 1);
                }

                // Compute the indirect irradiance, store it in delta_irradiance_texture and
                // accumulate it in irradiance_texture_.
                compute.SetTexture(compute_indirect_irradiance, "deltaIrradianceWrite", buffer.DeltaIrradianceTexture); //0
                compute.SetTexture(compute_indirect_irradiance, "irradianceWrite", buffer.IrradianceArray[WRITE]); //1
                compute.SetTexture(compute_indirect_irradiance, "irradianceRead", buffer.IrradianceArray[READ]);
                compute.SetTexture(compute_indirect_irradiance, "singleRayleighScatteringRead", buffer.DeltaRayleighScatteringTexture);
                compute.SetTexture(compute_indirect_irradiance, "singleMieScatteringRead", buffer.DeltaMieScatteringTexture);
                compute.SetTexture(compute_indirect_irradiance, "multipleScatteringRead", buffer.DeltaMultipleScatteringTexture);
                compute.SetInt("scatteringOrder", scattering_order - 1);
                compute.SetVector("blend", new Vector4(0, 1, 0, 0));

                compute.Dispatch(compute_indirect_irradiance, CONSTANTS.IRRADIANCE_WIDTH / NUM_THREADS, CONSTANTS.IRRADIANCE_HEIGHT / NUM_THREADS, 1);
                Swap(buffer.IrradianceArray);

                // Compute the multiple scattering, store it in
                // delta_multiple_scattering_texture, and accumulate it in
                // scattering_texture_.
                compute.SetTexture(compute_multiple_scattering, "deltaMultipleScatteringWrite", buffer.DeltaMultipleScatteringTexture); //0
                compute.SetTexture(compute_multiple_scattering, "scatteringWrite", buffer.ScatteringArray[WRITE]); //1
                compute.SetTexture(compute_multiple_scattering, "scatteringRead", buffer.ScatteringArray[READ]);
                compute.SetTexture(compute_multiple_scattering, "transmittanceRead", buffer.TransmittanceArray[READ]);
                compute.SetTexture(compute_multiple_scattering, "deltaScatteringDensityRead", buffer.DeltaScatteringDensityTexture);
                compute.SetVector("blend", new Vector4(0, 1, 0, 0));

                for (int layer = 0; layer < CONSTANTS.SCATTERING_DEPTH; ++layer) 
                {
                    compute.SetInt("layer", layer);
                    compute.Dispatch(compute_multiple_scattering, CONSTANTS.SCATTERING_WIDTH / NUM_THREADS, CONSTANTS.SCATTERING_HEIGHT / NUM_THREADS, 1);
                }
                Swap(buffer.ScatteringArray);

            }

            return;
        }

        private void Swap(RenderTexture[] arr)
        {
            RenderTexture tmp = arr[READ];
            arr[READ] = arr[WRITE];
            arr[WRITE] = tmp;
        }

        private void ReleaseTexture(RenderTexture tex)
        {
            if (tex == null) return;
            GameObject.DestroyImmediate(tex);
        }

        /// <summary>
        /// Convert a double array to a float matrix to bind to compute shader.
        /// Array is transposed.
        /// </summary>
        private float[] ToMatrix(double[] arr)
        {
            return new float[]
            {
                    (float)arr[0], (float)arr[3], (float)arr[6], 0,
                    (float)arr[1], (float)arr[4], (float)arr[7], 0,
                    (float)arr[2], (float)arr[5], (float)arr[8], 0,
                    0, 0, 0, 1
            };
        }

    }

}