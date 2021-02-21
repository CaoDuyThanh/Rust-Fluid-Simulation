mod fluid2d;
use fluid2d::Fluid2D;
use olc_pixel_game_engine as olc;

struct FluidSimulation {
    fluid2d: Fluid2D,
    sprite: olc::Sprite,

    pre_mouse: (i32, i32),
}

fn num_to_rgb(val: f32, max_val: i32) -> (u8, u8, u8) {
    let i = val / max_val as f32;
    let r = ((0.024 * i + 0.0).sin() * 127.0 + 128.0).round() as u8;
    let g = ((0.024 * i + 2.0).sin() * 127.0 + 128.0).round() as u8;
    let b = ((0.024 * i + 4.0).sin() * 127.0 + 128.0).round() as u8;
    (r,g,b)
}

impl olc::Application for FluidSimulation {
    fn on_user_create(&mut self) -> Result<(), olc::Error> {
        Ok(())
    }

    fn on_user_update(&mut self, _eslapse_time: f32) -> Result<(), olc::Error> {
        olc::clear(olc::BLACK);

        if olc::get_mouse(0).held {
            let (mx, my) = (olc::get_mouse_x(), olc::get_mouse_y());
            if self.pre_mouse.0 == -1 || self.pre_mouse.1 == -1 {
                self.pre_mouse = (mx, my);
            }
            else {
                self.fluid2d.add_density(mx as usize, my as usize, 2, 100.0);
                self.fluid2d.add_velocity(mx as usize, my as usize, 2, (mx - self.pre_mouse.0) as f32 * 1.0, (my - self.pre_mouse.1) as f32 * 1.0);
                self.pre_mouse = (mx, my);
            }
        }

        if olc::get_mouse(0).released {
            self.pre_mouse = (-1, -1);
        }

        self.fluid2d.step(1);

        // Update sprite
        for i in 0..self.fluid2d.size {
            for j in 0..self.fluid2d.size {
                let value = self.fluid2d.density[[i, j]];
                let (r, g, b) = num_to_rgb(value, 3);
                self.sprite.set_pixel(
                    i as i32,
                    j as i32,
                    olc::Pixel::rgb(r, g, b),
                    );
                // } else {
                //     self.sprite.set_pixel(
                //         i as i32,
                //         j as i32,
                //         olc::Pixel::rgba(255, 0, 255, value as u8),
                //     );
                // }
            }
        }

        olc::draw_sprite(0, 0, &self.sprite);
        // olc::draw_string(40, 40, "Hello, world!", olc::WHITE)?;

        Ok(())
    }

    fn on_user_destroy(&mut self) -> Result<(), olc::Error> {
        Ok(())
    }
}

fn main() {
    let mut example = FluidSimulation {
        fluid2d: Fluid2D::new(128, 0.001, 0.00001, 0.000001),
        sprite: olc::Sprite::with_dims(128, 128),
        pre_mouse: (-1, -1),
    };
    olc::start("Hello, World!", &mut example, 128, 128, 7, 7).unwrap();
}
