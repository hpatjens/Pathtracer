use common::*;

use scene::{Scene, Sky, Material, Sphere, Plane, PBRParameters};

#[derive(Clone, Debug, new)]
pub struct ParseError {
    message: String,
    position: TextPosition,
}

type Text<'a> = &'a str;

type ParseSuccess<'a, T> = (T, ParseContext<'a>);
type ParseResult<'a, T> = Result<ParseSuccess<'a, T>, ParseError>;

#[derive(Clone, Copy, Debug, new)]
pub struct TextPosition {
    line: usize,
    column: usize,
}

impl TextPosition {
    fn advanced_column(&mut self) {
        self.column += 1;
    }
    fn advanced_column_n(&mut self, n: usize) {
        self.column += n;
    }
    fn advanced_line(&mut self) {
        self.line += 1;
        self.column = 1;
    }
}

#[derive(Clone, Debug, new)]
struct ParseContext<'a> {
    text: Text<'a>,
    position: TextPosition,
}

fn success<T>(value: T, context: ParseContext) -> ParseResult<T> {
    Ok((value, context))
}

fn error<'a, T>(message: String, context_before_error: &ParseContext) -> ParseResult<'a, T> {
    Err(ParseError::new(message, context_before_error.position))
}

pub fn parse_scene(content: &str) -> Result<Scene, ParseError> {
    let context = ParseContext::new(content, TextPosition::new(1, 1));

    let mut running_context = context.clone();

    let mut spheres = Vec::new();
    let mut planes = Vec::new();
    let mut sky = Sky::Constant(Vec3::new(1.0, 1.0, 1.0));

    // @TODO: This structure is needed more often. There should be a function doing this.
    loop {
        let (_, context) = parse_whitespace(&running_context)?;
        if context.text.len() == 0 {
            break;
        }

        if let Ok((sphere, context)) = parse_sphere(&context) {
            spheres.push(sphere);
            running_context = context;
            continue;
        }

        if let Ok((plane, context)) = parse_plane(&context) {
            planes.push(plane);
            running_context = context;
            continue;
        }

        if let Ok((parsed_sky, context)) = parse_sky(&context) {
            sky = parsed_sky;
            running_context = context;
            continue;
        }

        return Err(ParseError::new(String::from("Expected sphere or plane."), context.position));
    }

    Ok(Scene::new(sky, spheres, planes))
}

enum SkyType { Constant, HDRI }

fn parse_whitespace_and_sky_type<'a>(context: &ParseContext<'a>) -> ParseResult<'a, SkyType> {
    let (_, context) = parse_whitespace(&context)?;

    if let Ok((_, context)) = parse_whitespace_and_string(&context, "constant") {
        return success(SkyType::Constant, context);
    }

    if let Ok((_, context)) = parse_whitespace_and_string(&context, "hdri") {
        return success(SkyType::HDRI, context);
    }

    error(String::from("Unknown sky type."), &context)
}

fn parse_whitespace_and_constant_sky<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Sky> {
    let (_       , context) = parse_whitespace_and_string(&context, "{")?;
    let (_       , context) = parse_whitespace_and_string(&context, "radiance")?;
    let (_       , context) = parse_whitespace_and_string(&context, "=")?;
    let (radiance, context) = parse_whitespace_and_vec3(&context)?;
    let (_       , context) = parse_whitespace_and_string(&context, "}")?;
    success(Sky::Constant(radiance), context)
}

fn parse_whitespace_and_hdri_sky<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Sky> {
    let (_   , context) = parse_whitespace_and_string(&context, "{")?;
    let (_   , context) = parse_whitespace_and_string(&context, "path")?;
    let (_   , context) = parse_whitespace_and_string(&context, "=")?;
    let (path, context) = parse_whitespace_and_path(&context)?;
    let (_   , context) = parse_whitespace_and_string(&context, "}")?;
    success(Sky::HDRI(path, None), context)
}

fn parse_sky<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Sky> {
    let (_, context) = parse_whitespace_and_string(&context, "sky")?;
    let (_, context) = parse_whitespace_and_string(&context, "{")?;

    let (sky_type, context) = parse_whitespace_and_sky_type(&context)?;

    let (sky, context) = match sky_type {
        SkyType::Constant => parse_whitespace_and_constant_sky(&context),
        SkyType::HDRI => parse_whitespace_and_hdri_sky(&context),
    }?;

    let (_, context) = parse_whitespace_and_string(&context, "}")?;
    
    success(sky, context)
}

fn parse_whitespace_and_vec3<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Vec3> {
    let (_, context) = parse_whitespace_and_string(&context, "(")?;

    // @TODO: There should be a function for parsing separated values of the same type.
    let (x, context) = parse_whitespace_and_f32(&context)?;
    let (_, context) = parse_whitespace_and_string(&context, ",")?;
    
    let (y, context) = parse_whitespace_and_f32(&context)?;
    let (_, context) = parse_whitespace_and_string(&context, ",")?;

    let (z, context) = parse_whitespace_and_f32(&context)?;
    let (_, context) = parse_whitespace_and_string(&context, ")")?;

    success(Vec3::new(x, y, z), context)
}

enum MaterialType { Physically, Emissive, Translucent }
fn parse_whitespace_and_material_type<'a>(context: &ParseContext<'a>) -> ParseResult<'a, MaterialType> {
    if let Ok((_, context)) = parse_whitespace_and_string(&context, "physically") {
        return success(MaterialType::Physically, context);
    }

    if let Ok((_, context)) = parse_whitespace_and_string(&context, "emissive") {
        return success(MaterialType::Emissive, context);
    }

    if let Ok((_, context)) = parse_whitespace_and_string(&context, "translucent") {
        return success(MaterialType::Translucent, context);
    }

    // @TODO: Add a default material for quickly setting up a scene.

    error(String::from("Unknown material type."), context)
}

fn parse_whitespace_and_physically_material<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Material> {
    let (_, context) = parse_whitespace_and_string(&context, "{")?;

    // @TODO: The order should not matter. This problem comes up more than once.
    let (_           , context) = parse_whitespace_and_string(&context, "reflectivity")?;
    let (_           , context) = parse_whitespace_and_string(&context, "=")?;
    let (reflectivity, context) = parse_whitespace_and_vec3(&context)?;

    let (_           , context) = parse_whitespace_and_newline(&context)?;

    let (_           , context) = parse_whitespace_and_string(&context, "roughness")?;
    let (_           , context) = parse_whitespace_and_string(&context, "=")?;
    let (roughness   , context) = parse_whitespace_and_f32(&context)?;

    let (_           , context) = parse_whitespace_and_newline(&context)?;

    let (_           , context) = parse_whitespace_and_string(&context, "metalness")?;
    let (_           , context) = parse_whitespace_and_string(&context, "=")?;
    let (metalness   , context) = parse_whitespace_and_f32(&context)?;

    let (_            , context) = parse_whitespace_and_string(&context, "}")?;

    let pbr_parameters = PBRParameters {
        reflectivity: reflectivity,
        roughness: roughness,
        metalness: metalness,
    };

    success(Material::Physically(pbr_parameters), context)
}

fn parse_whitespace_and_emissive_material<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Material> {
    let (_, context) = parse_whitespace_and_string(&context, "{")?;

    let (_       , context) = parse_whitespace_and_string(&context, "radiance")?;
    let (_       , context) = parse_whitespace_and_string(&context, "=")?;
    let (radiance, context) = parse_whitespace_and_vec3(&context)?;

    let (_       , context) = parse_whitespace_and_string(&context, "}")?;

    success(Material::Emissive(radiance), context)
}

fn parse_whitespace_and_translucent_material<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Material> {
    let (_, context) = parse_whitespace_and_string(&context, "{")?;

    let (_  , context) = parse_whitespace_and_string(&context, "ior")?;
    let (_  , context) = parse_whitespace_and_string(&context, "=")?;
    let (ior, context) = parse_whitespace_and_f32(&context)?;

    let (_  , context) = parse_whitespace_and_string(&context, "}")?;

    success(Material::Translucent(ior), context)
}

fn parse_whitespace_and_material<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Material> {
    let (material_type, context) = parse_whitespace_and_material_type(&context)?;

    match material_type {
        MaterialType::Physically => parse_whitespace_and_physically_material(&context),
        MaterialType::Emissive => parse_whitespace_and_emissive_material(&context),
        MaterialType::Translucent => parse_whitespace_and_translucent_material(&context),
    }
}

fn parse_sphere<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Sphere> {
    let (_, context) = parse_whitespace_and_string(&context, "sphere")?;
    let (_, context) = parse_whitespace_and_string(&context, "{")?;

    // @TODO: These closures are not needed.
    let origin = |context| {
        let (_, context) = parse_whitespace_and_string(&context, "origin")?;
        let (_, context) = parse_whitespace_and_string(&context, "=")?;
        let (v, context) = parse_whitespace_and_vec3(&context)?;
        success(v, context)
    };

    let radius = |context| {
        let (_, context) = parse_whitespace_and_string(&context, "radius")?;
        let (_, context) = parse_whitespace_and_string(&context, "=")?;
        let (r, context) = parse_whitespace_and_f32(&context)?;
        success(r, context)
    };

    let material = |context| {
        let (_, context) = parse_whitespace_and_string(&context, "material")?;
        let (_, context) = parse_whitespace_and_string(&context, "=")?;
        let (m, context) = parse_whitespace_and_material(&context)?;
        success(m, context)
    };

    // @TODO: The order should not matter. This problem comes up more than once.
    let (origin  , context) = origin(context)?;
    let (_       , context) = parse_whitespace_and_newline(&context)?;
    let (radius  , context) = radius(context)?;
    let (_       , context) = parse_whitespace_and_newline(&context)?;
    let (material, context) = material(context)?;
    let (_       , context) = parse_whitespace_and_string(&context, "}")?;

    success(Sphere::new(origin, radius, material), context)
}

fn parse_plane<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Plane> {
    let (_       , context) = parse_whitespace_and_string(&context, "plane")?;
    let (_       , context) = parse_whitespace_and_string(&context, "{")?;

    // @TODO: The order should not matter. This problem comes up more than once.
    let (_       , context) = parse_whitespace_and_string(&context, "origin")?;
    let (_       , context) = parse_whitespace_and_string(&context, "=")?;
    let (origin  , context) = parse_whitespace_and_vec3(&context)?;

    let (_       , context) = parse_whitespace_and_newline(&context)?;

    let (_       , context) = parse_whitespace_and_string(&context, "u")?;
    let (_       , context) = parse_whitespace_and_string(&context, "=")?;
    let (u       , context) = parse_whitespace_and_vec3(&context)?;

    let (_       , context) = parse_whitespace_and_newline(&context)?;

    let (_       , context) = parse_whitespace_and_string(&context, "v")?;
    let (_       , context) = parse_whitespace_and_string(&context, "=")?;
    let (v       , context) = parse_whitespace_and_vec3(&context)?;

    let (_       , context) = parse_whitespace_and_newline(&context)?;

    let (_       , context) = parse_whitespace_and_string(&context, "material")?;
    let (_       , context) = parse_whitespace_and_string(&context, "=")?;
    let (material, context) = parse_whitespace_and_material(&context)?;

    let (_       , context) = parse_whitespace_and_string(&context, "}")?;

    success(Plane::new(origin, u, v, material), context)
}

fn parse_whitespace_and_path<'a>(context: &ParseContext<'a>) -> ParseResult<'a, String> {
    let (_, context) = parse_whitespace(&context)?;

    // @TODO: How is a path defined?
    let ident = context.text.chars().take_while(|c| c.is_alphabetic() || c.is_numeric() || *c == '\\' || *c == '/' || *c == '_' || *c == '.');
    let count = ident.count();
    if count == 0 {
        error(String::from("Path expected."), &context)
    } else {
        let (path, rest) = context.text.split_at(count);
        let mut position = context.position;
        position.advanced_column_n(count);
        let new_context = ParseContext::new(rest, position);
        success(String::from(path), new_context)
    }
}

#[allow(dead_code)]
fn parse_whitespace_and_ident<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Text<'a>> {
    let (_    , context) = parse_whitespace(context)?;
    let (ident, context) = parse_alphabetic(&context)?;
    success(ident, context)
}

#[allow(dead_code)]
fn parse_alphabetic<'a>(context: &ParseContext<'a>) -> ParseResult<'a, Text<'a>> {
    let ident = context.text.chars().take_while(|c| c.is_alphabetic());
    let count = ident.count();
    if count == 0 {
        error(String::from("Identifier expected."), context)
    } else {
        let (ident, rest) = context.text.split_at(count);
        let mut position = context.position;
        position.advanced_column_n(count);
        let new_context = ParseContext::new(rest, position);
        success(ident, new_context)
    }
}

// Returns the number of newlines
fn parse_whitespace<'a>(context: &ParseContext<'a>) -> ParseResult<'a, usize> {
    let mut position = context.position;
    let count = {
        let whitespace = context.text.chars().take_while(|c| {
            if *c == '\n' {
                position.advanced_line();
                true
            } else if c.is_whitespace() {
                position.advanced_column();
                true
            } else {
                false
            }
        });
        whitespace.count()
    };
    let line_difference = position.line - context.position.line;
    let (_, rest) = context.text.split_at(count);
    let new_context = ParseContext::new(rest, position);
    success(line_difference, new_context)
}

fn parse_whitespace_and_string<'a>(context: &ParseContext<'a>, string: &str) -> ParseResult<'a, ()> {
    let (_, context) = parse_whitespace(context)?;
    let (_, context) = parse_string(&context, string)?;
    success((), context)
}

fn parse_string<'a>(context: &ParseContext<'a>, string: &str) -> ParseResult<'a, ()> {
    let starts_with = context.text.starts_with(string);
    if starts_with {
        let len = string.len();
        let (_, rest) = context.text.split_at(len);
        let mut position = context.position;
        position.advanced_column_n(len);
        let new_context = ParseContext::new(rest, position);
        success((), new_context)
    } else {
        error(format!("Could not parse the string \"{}\".", string), context)
    }
}

fn parse_sign<'a>(context: &ParseContext<'a>) -> ParseResult<'a, i8> {
    let negative = parse_string(context, "-");
    if let Some((_, context)) = negative.ok() {
        return success(-1, context);
    }

    let positive = parse_string(context, "+");
    if let Some((_, context)) = positive.ok() {
        return success(1, context);
    }

    // Note that this is the original context (not the same as the return above)
    success(1, context.clone())
}

fn parse_u32<'a>(context: &ParseContext<'a>) -> ParseResult<'a, u32> {
    let numbers = context.text.chars().take_while(|c| c.is_numeric());
    let count = numbers.count();
    if count > 0 {
        let (numbers, rest) = context.text.split_at(count);
        let value = numbers.parse().unwrap();
        let mut position = context.position;
        position.advanced_column_n(count);
        let new_context = ParseContext::new(rest, position);
        success(value, new_context)
    } else {
        error(String::from("Could not parse the u32."), &context)
    }
}

#[allow(dead_code)]
fn parse_whitespace_and_i32<'a>(context: &ParseContext<'a>) -> ParseResult<'a, i32> {
    let (_    , context) = parse_whitespace(context)?;
    let (value, context) = parse_i32(&context)?;
    success(value, context)
}

fn parse_i32<'a>(context: &ParseContext<'a>) -> ParseResult<'a, i32> {
    let (sign, context) = parse_sign(context)?;

    let numbers = context.text.chars().take_while(|c| c.is_numeric());
    let count = numbers.count();
    if count > 0 {
        let (numbers, rest) = context.text.split_at(count);
        let value = numbers.parse::<i32>().unwrap();
        let mut position = context.position;
        position.advanced_column_n(count);
        let new_context = ParseContext::new(rest, position);
        success(sign as i32*value, new_context)
    } else {
        error(String::from("Could not parse the i32."), &context)
    }
}

fn parse_f32<'a>(context: &ParseContext<'a>) -> ParseResult<'a, f32> {
    let ref start_context = context;

    let (_, context) = parse_sign(context)?;
    let (_, context) = parse_i32(&context)?;
    let (_, context) = parse_string(&context, ".")?;
    let (_, context) = parse_u32(&context)?;

    let len = context.position.column - start_context.position.column;

    let (characters, rest) = start_context.text.split_at(len);
    let value = characters.parse().unwrap();

    let mut position = start_context.position;
    position.advanced_column_n(len);
    let new_context = ParseContext::new(rest, position);

    success(value, new_context)
}

fn parse_whitespace_and_f32<'a>(context: &ParseContext<'a>) -> ParseResult<'a, f32> {
    let (_    , context) = parse_whitespace(context)?;
    let (value, context) = parse_f32(&context)?;
    success(value, context)
}

fn parse_whitespace_and_newline<'a>(context: &ParseContext<'a>) -> ParseResult<'a, ()> {
    let (num_newlines, new_context) = parse_whitespace(&context)?;
    if num_newlines == 0 {
        error(String::from("Expected at least one newline."), context)
    } else {
        success((), new_context)
    }
}